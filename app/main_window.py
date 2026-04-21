import json
import os
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Dict, List, Optional

from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from app.services.processor import process_mutagenesis


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SETTINGS_PATH = PROJECT_ROOT / "resources" / "settings.json"
RESULT_COLUMNS = [
    "variant pdb",
    "sequence",
    "EvoEF2 binding",
    "Interface area",
    "Interface solvation E",
    "p-value",
    "status",
]
INPUT_RESULT_FIELDS = [
    ("input pdb", "Input PDB"),
    ("residue range", "Residues"),
    ("sequence", "Sequence"),
    ("EvoEF2 binding", "EvoEF2"),
    ("Interface area", "PISA Area"),
    ("Interface solvation E", "PISA Solv E"),
    ("p-value", "PISA p-value"),
]
MUTAGENESIS_RESULT_COLUMNS = [
    "mutation",
    "count",
    "EvoEF2 improved",
    "Interface area improved",
    "Solvation E improved",
]


def default_settings() -> Dict[str, str]:
    return {
        "proteinmpnn_dir": "",
        "proteinmpnn_env_setup": "",
        "proteinmpnn_seed": "37",
        "proteinmpnn_batch_size": "1",
        "evoef2_executable": "EvoEF2",
        "evoef2_threads": "1",
    }


def load_settings() -> Dict[str, str]:
    settings = default_settings()
    if SETTINGS_PATH.exists():
        try:
            loaded = json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                for key in settings:
                    settings[key] = str(loaded.get(key, settings[key])).strip()
        except Exception:
            pass
    return settings


def save_settings(settings: Dict[str, str]) -> None:
    SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(settings, indent=2), encoding="utf-8")


class SettingsDialog(QDialog):
    def __init__(self, settings: Dict[str, str], parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.resize(820, 420)
        self._build_ui(settings)

    def _build_ui(self, settings: Dict[str, str]) -> None:
        layout = QVBoxLayout(self)
        form = QFormLayout()
        form.setSpacing(10)

        self.proteinmpnn_dir_edit = self._path_row(
            form,
            "ProteinMPNN directory",
            settings.get("proteinmpnn_dir", ""),
            self._browse_proteinmpnn_dir,
            directory=True,
        )
        self.proteinmpnn_env_edit = QPlainTextEdit(settings.get("proteinmpnn_env_setup", ""))
        self.proteinmpnn_env_edit.setPlaceholderText("Optional shell setup, for example module load and conda activate commands")
        self.proteinmpnn_env_edit.setFixedHeight(86)
        form.addRow("ProteinMPNN env setup", self.proteinmpnn_env_edit)
        self.proteinmpnn_seed_edit = QLineEdit(settings.get("proteinmpnn_seed", "37"))
        form.addRow("ProteinMPNN seed", self.proteinmpnn_seed_edit)
        self.proteinmpnn_batch_edit = QLineEdit(settings.get("proteinmpnn_batch_size", "1"))
        form.addRow("ProteinMPNN batch size", self.proteinmpnn_batch_edit)
        self.evoef2_edit = self._path_row(
            form,
            "EvoEF2 executable",
            settings.get("evoef2_executable", "EvoEF2"),
            self._browse_evoef2,
            directory=False,
        )
        self.evoef2_threads_edit = QLineEdit(settings.get("evoef2_threads", "1"))
        form.addRow("EvoEF2 threads", self.evoef2_threads_edit)

        layout.addLayout(form)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _path_row(self, form: QFormLayout, label: str, value: str, slot, directory: bool) -> QLineEdit:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        edit = QLineEdit(value)
        button = QPushButton("Browse")
        button.clicked.connect(slot)
        layout.addWidget(edit, 1)
        layout.addWidget(button)
        form.addRow(label, row)
        return edit

    def _browse_proteinmpnn_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select ProteinMPNN directory", self.proteinmpnn_dir_edit.text())
        if path:
            self.proteinmpnn_dir_edit.setText(path)

    def _browse_evoef2(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select EvoEF2 executable", self.evoef2_edit.text())
        if path:
            self.evoef2_edit.setText(path)

    def values(self) -> Dict[str, str]:
        return {
            "proteinmpnn_dir": self.proteinmpnn_dir_edit.text().strip(),
            "proteinmpnn_env_setup": self.proteinmpnn_env_edit.toPlainText().strip(),
            "proteinmpnn_seed": self.proteinmpnn_seed_edit.text().strip() or "37",
            "proteinmpnn_batch_size": self.proteinmpnn_batch_edit.text().strip() or "1",
            "evoef2_executable": self.evoef2_edit.text().strip() or "EvoEF2",
            "evoef2_threads": self.evoef2_threads_edit.text().strip() or "1",
        }


class WorkerSignals(QObject):
    log = Signal(str)
    input_result = Signal(dict)
    result = Signal(dict)
    mutagenesis_result = Signal(list)
    progress = Signal(int)
    status = Signal(str)
    error = Signal(str)
    finished = Signal()


class MutagenesisWorker(QObject):
    def __init__(
        self,
        input_pdb: str,
        output_dir: str,
        target_chain: str,
        design_chain: str,
        residue_range: str,
        replicates: int,
        settings: Dict[str, str],
    ) -> None:
        super().__init__()
        self.input_pdb = input_pdb
        self.output_dir = output_dir
        self.target_chain = target_chain
        self.design_chain = design_chain
        self.residue_range = residue_range
        self.replicates = replicates
        self.settings = settings
        self.signals = WorkerSignals()
        self._cancel_requested = False
        self._processes: List[subprocess.Popen] = []

    def request_cancel(self) -> None:
        self._cancel_requested = True
        for process in list(self._processes):
            if process.poll() is None:
                try:
                    process.terminate()
                except Exception:
                    pass

    def run(self) -> None:
        try:
            self.signals.status.emit("Running")
            process_mutagenesis(
                input_pdb=self.input_pdb,
                output_dir=self.output_dir,
                target_chain=self.target_chain,
                design_chain=self.design_chain,
                residue_range=self.residue_range,
                replicates=self.replicates,
                settings=self.settings,
                log_callback=self.signals.log.emit,
                input_result_callback=self.signals.input_result.emit,
                result_callback=self.signals.result.emit,
                mutagenesis_result_callback=self.signals.mutagenesis_result.emit,
                progress_callback=self.signals.progress.emit,
                should_cancel_callback=lambda: self._cancel_requested,
                register_process_callback=self._register_process,
                unregister_process_callback=self._unregister_process,
            )
            self.signals.status.emit("Completed")
        except Exception as exc:
            if str(exc) == "Cancelled":
                self.signals.status.emit("Cancelled")
            else:
                self.signals.error.emit("".join(traceback.format_exception_only(type(exc), exc)).strip())
                self.signals.log.emit(traceback.format_exc())
                self.signals.status.emit("Failed")
        finally:
            self.signals.finished.emit()

    def _register_process(self, process: subprocess.Popen) -> None:
        self._processes.append(process)

    def _unregister_process(self, process: subprocess.Popen) -> None:
        if process in self._processes:
            self._processes.remove(process)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Binder Mutagenesis GUI")
        self.resize(1320, 820)
        self.settings = load_settings()
        self.worker: Optional[MutagenesisWorker] = None
        self.thread: Optional[QThread] = None
        self.rows_by_key: Dict[str, int] = {}
        self.input_result_values: Dict[str, QLabel] = {}
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(14, 14, 14, 14)
        root.setSpacing(10)

        title_row = QHBoxLayout()
        title = QLabel("Binder Mutagenesis GUI")
        title.setStyleSheet("font-size: 22px; font-weight: 700;")
        title_row.addWidget(title)
        title_row.addStretch(1)
        self.settings_btn = QPushButton("Settings")
        self.settings_btn.setMinimumHeight(34)
        self.settings_btn.clicked.connect(self._open_settings)
        title_row.addWidget(self.settings_btn)
        root.addLayout(title_row)

        top_row = QHBoxLayout()
        top_row.setSpacing(10)
        top_row.addWidget(self._build_input_group(), 3)
        top_row.addWidget(self._build_run_group(), 1)
        root.addLayout(top_row)

        root.addWidget(self._build_input_result_group())

        result_row = QHBoxLayout()
        result_row.setSpacing(10)
        result_row.addWidget(self._build_replication_output_group(), 3)
        result_row.addWidget(self._build_mutagenesis_result_group(), 2)
        root.addLayout(result_row, 1)

        log_group = QGroupBox("Log")
        log_layout = QVBoxLayout(log_group)
        self.log_edit = QPlainTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumBlockCount(4000)
        self.log_edit.setMaximumHeight(180)
        log_layout.addWidget(self.log_edit)
        root.addWidget(log_group)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setFormat("%p%")
        root.addWidget(self.progress)

        self.setCentralWidget(central)

    def _build_input_group(self) -> QGroupBox:
        input_group = QGroupBox("Input")
        input_grid = QGridLayout(input_group)
        input_grid.setContentsMargins(10, 12, 10, 10)
        input_grid.setHorizontalSpacing(8)
        input_grid.setVerticalSpacing(8)
        self.input_pdb_edit = QLineEdit()
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Optional. Leave blank to save in the input PDB directory")
        self.residue_range_edit = QLineEdit()
        self.residue_range_edit.setPlaceholderText("Example: 45-80 or 45,48,52-60")
        self.replicates_spin = QSpinBox()
        self.replicates_spin.setRange(1, 10000)
        self.replicates_spin.setValue(5000)

        input_grid.addWidget(QLabel("Input PDB"), 0, 0)
        input_grid.addWidget(self.input_pdb_edit, 0, 1)
        browse_pdb = QPushButton("Browse")
        browse_pdb.clicked.connect(self._browse_input_pdb)
        input_grid.addWidget(browse_pdb, 0, 2)
        input_grid.addWidget(QLabel("Output directory"), 1, 0)
        input_grid.addWidget(self.output_dir_edit, 1, 1)
        browse_output = QPushButton("Browse")
        browse_output.setFixedWidth(80)
        browse_output.clicked.connect(self._browse_output_dir)
        input_grid.addWidget(browse_output, 1, 2)

        options = ["", "A", "B", "C", "D", "E", "Other"]
        input_grid.addWidget(QLabel("Target"), 2, 0)
        target_row = QHBoxLayout()
        target_row.setContentsMargins(0, 0, 0, 0)
        target_row.setSpacing(6)
        self.target_chain_combo = QComboBox()
        self.target_chain_combo.addItems(options)
        self.target_chain_combo.setFixedWidth(120)
        self.target_chain_combo.setCurrentIndex(0)
        self.target_chain_other_edit = QLineEdit()
        self.target_chain_other_edit.setPlaceholderText("Enter target chain")
        self.target_chain_other_edit.setFixedWidth(140)
        self.target_chain_other_edit.setVisible(False)
        target_row.addWidget(self.target_chain_combo)
        target_row.addWidget(self.target_chain_other_edit)
        target_row.addStretch(1)
        target_row_widget = QWidget()
        target_row_widget.setLayout(target_row)
        input_grid.addWidget(target_row_widget, 2, 1, 1, 2)

        input_grid.addWidget(QLabel("Design"), 3, 0)
        design_row = QHBoxLayout()
        design_row.setContentsMargins(0, 0, 0, 0)
        design_row.setSpacing(6)
        self.design_chain_combo = QComboBox()
        self.design_chain_combo.addItems(options)
        self.design_chain_combo.setFixedWidth(120)
        self.design_chain_combo.setCurrentIndex(0)
        self.design_chain_other_edit = QLineEdit()
        self.design_chain_other_edit.setPlaceholderText("Enter design chain")
        self.design_chain_other_edit.setFixedWidth(140)
        self.design_chain_other_edit.setVisible(False)
        design_row.addWidget(self.design_chain_combo)
        design_row.addWidget(self.design_chain_other_edit)
        design_row.addStretch(1)
        design_row_widget = QWidget()
        design_row_widget.setLayout(design_row)
        input_grid.addWidget(design_row_widget, 3, 1, 1, 2)

        self.target_chain_combo.currentTextChanged.connect(self._on_target_chain_changed)
        self.design_chain_combo.currentTextChanged.connect(self._on_design_chain_changed)
        input_grid.addWidget(QLabel("Residues range to screen"), 4, 0)
        input_grid.addWidget(self.residue_range_edit, 4, 1)
        input_grid.addWidget(QLabel("ProteinMPNN replicates"), 5, 0)
        input_grid.addWidget(self.replicates_spin, 5, 1)
        return input_group

    def _build_run_group(self) -> QGroupBox:
        group = QGroupBox("Parameters")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(10, 12, 10, 10)
        layout.setSpacing(8)
        form = QFormLayout()
        form.setSpacing(6)
        self.proteinmpnn_seed_value = QLabel("")
        self.proteinmpnn_batch_value = QLabel("")
        self.evoef2_threads_value = QLabel("")
        for value_label in (
            self.proteinmpnn_seed_value,
            self.proteinmpnn_batch_value,
            self.evoef2_threads_value,
        ):
            value_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            value_label.setStyleSheet("color: #101828; font-weight: 600;")
        form.addRow("ProteinMPNN seed", self.proteinmpnn_seed_value)
        form.addRow("ProteinMPNN batch size", self.proteinmpnn_batch_value)
        form.addRow("EvoEF2 threads", self.evoef2_threads_value)
        layout.addLayout(form)
        layout.addStretch(1)
        self.run_btn = QPushButton("Run")
        self.run_btn.setMinimumWidth(180)
        self.run_btn.setMinimumHeight(44)
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #d92d20;
                color: white;
                font-size: 20px;
                font-weight: 700;
                border: 1px solid #b42318;
                border-radius: 8px;
                padding: 8px 18px;
            }
            QPushButton:hover {
                background-color: #b42318;
            }
            QPushButton:pressed {
                background-color: #912018;
            }
            QPushButton:disabled {
                background-color: #f1a7a1;
                color: #fff5f4;
                border-color: #f1a7a1;
            }
            """
        )
        self.run_btn.clicked.connect(self._start_run)
        layout.addWidget(self.run_btn, 0, Qt.AlignCenter)
        layout.addStretch(1)
        self._update_parameter_summary()
        return group

    def _update_parameter_summary(self) -> None:
        self.proteinmpnn_seed_value.setText(str(self.settings.get("proteinmpnn_seed", "")))
        self.proteinmpnn_batch_value.setText(str(self.settings.get("proteinmpnn_batch_size", "")))
        self.evoef2_threads_value.setText(str(self.settings.get("evoef2_threads", "")))

    def _build_input_result_group(self) -> QGroupBox:
        group = QGroupBox("Input Result")
        layout = QGridLayout(group)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(3)
        columns_per_row = 3
        for index, (key, label_text) in enumerate(INPUT_RESULT_FIELDS):
            row = index // columns_per_row
            column = index % columns_per_row
            field = QWidget()
            field_layout = QHBoxLayout(field)
            field_layout.setContentsMargins(0, 0, 0, 0)
            field_layout.setSpacing(4)
            label = QLabel(label_text)
            label.setMinimumWidth(86)
            label.setStyleSheet("color: #475467; font-weight: 700;")
            value = QLabel("")
            value.setTextInteractionFlags(Qt.TextSelectableByMouse)
            value.setStyleSheet("color: #101828;")
            value.setMinimumWidth(90)
            value.setWordWrap(False)
            field_layout.addWidget(label)
            field_layout.addWidget(value, 1)
            layout.addWidget(field, row, column)
            layout.setColumnStretch(column, 1)
            self.input_result_values[key] = value
        return group

    def _build_mutagenesis_result_group(self) -> QGroupBox:
        group = QGroupBox("Mutagenesis Result")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        self.mutagenesis_table = self._create_mutagenesis_table()
        layout.addWidget(self.mutagenesis_table)
        return group

    def _build_replication_output_group(self) -> QGroupBox:
        group = QGroupBox("Replication Output")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        self.table = QTableWidget(0, len(RESULT_COLUMNS))
        self.table.setHorizontalHeaderLabels(RESULT_COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        layout.addWidget(self.table)
        return group

    def _create_mutagenesis_table(self) -> QTableWidget:
        table = QTableWidget(0, len(MUTAGENESIS_RESULT_COLUMNS))
        table.setHorizontalHeaderLabels(MUTAGENESIS_RESULT_COLUMNS)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        table.horizontalHeader().setStretchLastSection(True)
        table.verticalHeader().setVisible(False)
        table.setAlternatingRowColors(True)
        table.setMinimumHeight(110)
        return table

    def _on_target_chain_changed(self, value: str) -> None:
        self.target_chain_other_edit.setVisible(value == "Other")

    def _on_design_chain_changed(self, value: str) -> None:
        self.design_chain_other_edit.setVisible(value == "Other")

    def _get_target_chain(self) -> str:
        value = self.target_chain_combo.currentText().strip()
        if value == "Other":
            return self.target_chain_other_edit.text().strip()
        return value

    def _get_design_chain(self) -> str:
        value = self.design_chain_combo.currentText().strip()
        if value == "Other":
            return self.design_chain_other_edit.text().strip()
        return value

    def _browse_input_pdb(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select input PDB", self.input_pdb_edit.text(), "PDB files (*.pdb)")
        if path:
            self.input_pdb_edit.setText(path)

    def _browse_output_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output directory", self.output_dir_edit.text())
        if path:
            self.output_dir_edit.setText(path)

    def _open_settings(self) -> None:
        dialog = SettingsDialog(self.settings, self)
        if dialog.exec() == QDialog.Accepted:
            self.settings = dialog.values()
            save_settings(self.settings)
            self._update_parameter_summary()

    def _start_run(self) -> None:
        try:
            input_pdb = self.input_pdb_edit.text().strip()
            output_dir = self.output_dir_edit.text().strip()
            target_chain = self._get_target_chain()
            design_chain = self._get_design_chain()
            residue_range = self.residue_range_edit.text().strip()
            replicates = self.replicates_spin.value()
            if not input_pdb:
                raise ValueError("Input PDB is required.")
            if not output_dir:
                output_dir = str(Path(input_pdb).expanduser().resolve().parent / "_binder_mutagenesis_gui")
            if not residue_range:
                raise ValueError("Residue range is required.")
        except Exception as exc:
            QMessageBox.warning(self, "Input error", str(exc))
            return

        self.table.setRowCount(0)
        self.mutagenesis_table.setRowCount(0)
        self.rows_by_key.clear()
        self._clear_input_result()
        self.log_edit.clear()
        self.progress.setValue(0)
        self.run_btn.setEnabled(False)
        self.settings_btn.setEnabled(False)

        self.worker = MutagenesisWorker(
            input_pdb=input_pdb,
            output_dir=output_dir,
            target_chain=target_chain,
            design_chain=design_chain,
            residue_range=residue_range,
            replicates=replicates,
            settings=dict(self.settings),
        )
        self.thread = QThread(self)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.signals.log.connect(self._append_log)
        self.worker.signals.input_result.connect(self._update_input_result)
        self.worker.signals.result.connect(self._upsert_result)
        self.worker.signals.mutagenesis_result.connect(self._replace_mutagenesis_result)
        self.worker.signals.progress.connect(self.progress.setValue)
        self.worker.signals.status.connect(self._append_status)
        self.worker.signals.error.connect(self._show_error)
        self.worker.signals.finished.connect(self._run_finished)
        self.worker.signals.finished.connect(self.thread.quit)
        self.worker.signals.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self._thread_finished)
        self.thread.start()

    def _append_log(self, text: str) -> None:
        self.log_edit.appendPlainText(text)

    def _append_status(self, text: str) -> None:
        if text:
            self._append_log(f"[STATUS] {text}")

    def _show_error(self, text: str) -> None:
        QMessageBox.critical(self, "Run failed", text)

    def _clear_input_result(self) -> None:
        for value_label in self.input_result_values.values():
            value_label.setText("")

    def _update_input_result(self, row: Dict[str, object]) -> None:
        for key, value_label in self.input_result_values.items():
            value = row.get(key, "")
            value_label.setText(str(value))
            value_label.setToolTip(str(value))

    def _upsert_result(self, row: Dict[str, object]) -> None:
        key = str(row.get("variant pdb", ""))
        if key in self.rows_by_key:
            table_row = self.rows_by_key[key]
        else:
            table_row = self.table.rowCount()
            self.table.insertRow(table_row)
            self.rows_by_key[key] = table_row
        for column_index, column_name in enumerate(RESULT_COLUMNS):
            item = QTableWidgetItem(str(row.get(column_name, "")))
            self.table.setItem(table_row, column_index, item)

    def _replace_mutagenesis_result(self, rows: List[Dict[str, object]]) -> None:
        self.mutagenesis_table.setRowCount(0)
        self._fill_mutagenesis_table(self.mutagenesis_table, rows)

    def _fill_mutagenesis_table(self, table: QTableWidget, rows: List[Dict[str, object]]) -> None:
        for row in rows:
            table_row = table.rowCount()
            table.insertRow(table_row)
            for column_index, column_name in enumerate(MUTAGENESIS_RESULT_COLUMNS):
                item = QTableWidgetItem(str(row.get(column_name, "")))
                if column_name != "mutation":
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                table.setItem(table_row, column_index, item)

    def _run_finished(self) -> None:
        self.run_btn.setEnabled(True)
        self.settings_btn.setEnabled(True)
        self.worker = None

    def _thread_finished(self) -> None:
        self.thread = None

    def closeEvent(self, event) -> None:
        if self.worker is not None and self.run_btn is not None and not self.run_btn.isEnabled():
            response = QMessageBox.question(
                self,
                "Exit",
                "A job is still running. Close the GUI and terminate active jobs?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if response != QMessageBox.Yes:
                event.ignore()
                return
            self.worker.request_cancel()
        event.accept()


def main() -> None:
    os.environ.setdefault("NO_AT_BRIDGE", "1")
    os.environ.setdefault("QT_LOGGING_RULES", "qt.qpa.theme.dbus=false;qt.qpa.theme.gnome=false")
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
