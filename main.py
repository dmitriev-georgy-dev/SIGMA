"""Sequence Interface for Genetic Modeling and Analysis (SIGMA) v1.2"""

import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QTextEdit, QVBoxLayout,
    QHBoxLayout, QWidget, QPushButton, QSplitter, QComboBox, QLabel
)
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtCore import Qt, QThread

# Модули SIGMA v1.2
from pdb_finder import LocalPDBFinder
from protein_db import PROTEIN_REGISTRY
from seq_analyzer import SequenceAnalyzer, VALID_DNA_CHARS, VALID_PROTEIN_CHARS
from syntax_tools import SequenceSyntaxHighlighter, highlight_mutations, clear_extra_highlights
from viewer_3d import Protein3DViewer

# Принудительно отключаем черные списки видеокарт Chromium на уровне ОС
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--ignore-gpu-blocklist --disable-gpu-driver-bug-workarounds"
os.environ["QT_ANGLE_PLATFORM"] = "d3d11"

FREE_INPUT_LABEL = "— Свободный ввод (без эталона) —"


def get_resource_path(relative_path):
    """Динамически определяет путь к ресурсам (работает и в .py, и в собранном .exe)"""
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


if not hasattr(sys, '_MEIPASS'):
    if not os.path.exists('structures'):
        os.makedirs('structures')


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SIGMA")
        self.resize(1180, 700)

        # Инициализируем движки анализа и 3D-визуализации
        self.analyzer = SequenceAnalyzer()
        self.viewer_3d = Protein3DViewer()
        self.current_protein_mutations = []

        self.setup_ui()

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        window_layout = QVBoxLayout(main_widget)
        window_layout.setContentsMargins(6, 6, 6, 6)

        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # --- ЛЕВАЯ ЧАСТЬ (Редакторы ДНК и Белка) ---
        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        left_layout.setContentsMargins(0, 0, 4, 0)

        # Верхняя панель выбора белка + кнопка сброса
        selector_row = QHBoxLayout()
        self.protein_selector = QComboBox()
        self.protein_selector.addItems(PROTEIN_REGISTRY.keys())
        self.protein_selector.addItem(FREE_INPUT_LABEL)
        self.protein_selector.currentIndexChanged.connect(self.on_protein_selected)

        self.reset_button = QPushButton("Сбросить")
        self.reset_button.setFixedWidth(90)
        self.reset_button.setToolTip("Вернуть эталонные последовательности ДНК и белка")
        self.reset_button.clicked.connect(self.on_protein_selected)

        selector_row.addWidget(self.protein_selector, stretch=1)
        selector_row.addWidget(self.reset_button)

        left_layout.addWidget(QLabel("Выберите эталонный белок или режим свободного ввода:"))
        left_layout.addLayout(selector_row)

        # Разделитель для двух уровней редактирования
        text_splitter = QSplitter(Qt.Orientation.Vertical)

        # Блок 1: Низкий уровень (ДНК)
        dna_panel = QWidget()
        dna_layout = QVBoxLayout(dna_panel)
        dna_layout.setContentsMargins(0, 4, 0, 4)

        self.dna_editor = QTextEdit()
        self.dna_editor.setPlaceholderText("Введите кодирующую ДНК (только символы A, T, G, C)...")
        self.dna_highlighter = SequenceSyntaxHighlighter(self.dna_editor.document(), VALID_DNA_CHARS)

        self.dna_run_button = QPushButton("Транслировать ДНК и смоделировать 3D")
        self.dna_run_button.clicked.connect(self.run_from_dna)

        dna_layout.addWidget(QLabel("Кодирующая цепь ДНК (низкий уровень):"))
        dna_layout.addWidget(self.dna_editor)
        dna_layout.addWidget(self.dna_run_button)

        # Блок 2: Высокий уровень (Белок)
        protein_panel = QWidget()
        protein_layout = QVBoxLayout(protein_panel)
        protein_layout.setContentsMargins(0, 4, 0, 0)

        self.protein_reader = QTextEdit()
        self.protein_reader.setPlaceholderText("Аминокислотная последовательность (20 стандартных аминокислот)...")
        self.protein_highlighter = SequenceSyntaxHighlighter(self.protein_reader.document(), VALID_PROTEIN_CHARS)

        self.protein_run_button = QPushButton("Анализировать белок и смоделировать 3D")
        self.protein_run_button.clicked.connect(self.run_from_protein)

        protein_layout.addWidget(QLabel("Аминокислотная цепь белка (высокий уровень):"))
        protein_layout.addWidget(self.protein_reader)
        protein_layout.addWidget(self.protein_run_button)

        text_splitter.addWidget(dna_panel)
        text_splitter.addWidget(protein_panel)
        left_layout.addWidget(text_splitter)

        # --- ПРАВАЯ ЧАСТЬ (Интерактивная 3D Графика) ---
        self.graphics_panel = QWebEngineView()
        self.set_placeholder_3d("Ожидание запуска анализа и 3D-моделирования...")

        main_splitter.addWidget(left_container)
        main_splitter.addWidget(self.graphics_panel)
        main_splitter.setSizes([560, 620])

        # --- НИЖНЯЯ ЧАСТЬ (Строка логов) ---
        self.status_log = QLabel("[Система] SIGMA v1.2 готова к работе.")
        self.status_log.setStyleSheet("color: #ddd; background-color: #252525; padding: 6px; font-family: monospace;")
        self.status_log.setFixedHeight(32)  # <-- Фиксируем аккуратную высоту строки

        window_layout.addWidget(main_splitter, stretch=1)  # <-- Отдаем всё пространство рабочей области
        window_layout.addWidget(self.status_log, stretch=0)

        # Загружаем первый белок при старте
        self.on_protein_selected()

    def set_placeholder_3d(self, message: str):
        self.graphics_panel.setHtml(
            f"<body style='background-color:#262626; color:#888; display:flex; justify-content:center; "
            f"align-items:center; height:100vh; margin:0; font-family:sans-serif;'>{message}</body>"
        )

    def _get_current_references(self) -> tuple[str, str, str]:
        """Возвращает (имя_белка, эталонная_ДНК, эталонный_белок) для выбранного пункта списка"""
        selected = self.protein_selector.currentText()
        if selected in PROTEIN_REGISTRY:
            data = PROTEIN_REGISTRY[selected]
            return selected, data["dna"], data["protein"]
        return "", "", ""

    def on_protein_selected(self):
        """Подставляет эталонные ДНК и белок или очищает поля в свободном режиме"""
        selected_name = self.protein_selector.currentText()
        self.current_protein_mutations = []

        self.dna_editor.blockSignals(True)
        self.protein_reader.blockSignals(True)

        if selected_name in PROTEIN_REGISTRY:
            data = PROTEIN_REGISTRY[selected_name]
            self.dna_editor.setPlainText(data["dna"])
            self.protein_reader.setPlainText(data["protein"])
            clear_extra_highlights(self.dna_editor)
            clear_extra_highlights(self.protein_reader)
            self.status_log.setText(f"[Система] Загружен эталон: {selected_name}")
        else:
            self.dna_editor.clear()
            self.protein_reader.clear()
            self.status_log.setText("[Система] Режим свободного ввода. Вставьте ДНК или белок для поиска в базе.")

        self.dna_editor.blockSignals(False)
        self.protein_reader.blockSignals(False)

    def run_from_dna(self):
        """Обработка нажатия кнопки под окном ДНК"""
        base_name, ref_dna, ref_protein = self._get_current_references()
        raw_dna = self.dna_editor.toPlainText()

        result = self.analyzer.analyze_from_dna(raw_dna, ref_dna, ref_protein)

        if not result.is_valid:
            self.status_log.setText(f"[Ошибка валидации] {result.error_message}")
            return

        # Обновляем нижнее поле с белком и подсвечиваем точечные мутации в обоих полях
        self.protein_reader.blockSignals(True)
        self.protein_reader.setPlainText(result.protein_seq)
        self.protein_reader.blockSignals(False)

        highlight_mutations(self.dna_editor, result.dna_mutations, is_protein=False)
        highlight_mutations(self.protein_reader, result.protein_mutations, is_protein=True)

        self.current_protein_mutations = result.protein_mutations
        self.status_log.setText(f"[Анализатор ДНК] {result.summary_log}")

        # Проверяем, разрешено ли строить 3D-модель
        if not result.can_render_3d:
            self.graphics_panel.setHtml(self.viewer_3d.build_blocked_html(result.render_block_reason))
            return

        self._start_pdb_search(result.protein_seq, base_name)

    def run_from_protein(self):
        """Обработка нажатия кнопки под окном Белка"""
        base_name, _, ref_protein = self._get_current_references()
        raw_protein = self.protein_reader.toPlainText()

        result = self.analyzer.analyze_from_protein(raw_protein, ref_protein)

        if not result.is_valid:
            self.status_log.setText(f"[Ошибка валидации] {result.error_message}")
            return

        # При прямой правке белка сбрасываем подсветку ДНК (так как генетический код вырожден)
        clear_extra_highlights(self.dna_editor)
        highlight_mutations(self.protein_reader, result.protein_mutations, is_protein=True)

        self.current_protein_mutations = result.protein_mutations
        self.status_log.setText(f"[Анализатор белка] {result.summary_log}")

        if not result.can_render_3d:
            self.graphics_panel.setHtml(self.viewer_3d.build_blocked_html(result.render_block_reason))
            return

        self._start_pdb_search(result.protein_seq, base_name)

    def _start_pdb_search(self, protein_seq: str, base_name: str):
        """Запускает фоновый поиск файла структуры в отдельном потоке"""
        self.thread = QThread()
        self.finder = LocalPDBFinder(protein_seq, get_resource_path, base_name)
        self.finder.moveToThread(self.thread)

        self.thread.started.connect(self.finder.run)
        self.finder.finished.connect(self.on_search_success)
        self.finder.error.connect(self.on_search_error)
        self.finder.finished.connect(self.thread.quit)
        self.finder.error.connect(self.thread.quit)
        self.finder.finished.connect(self.finder.deleteLater)
        self.finder.error.connect(self.finder.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()

    def on_search_success(self, file_path: str):
        try:
            html_content = self.viewer_3d.build_protein_html(file_path, self.current_protein_mutations)

            self.graphics_panel.setHtml(html_content)
        except Exception as e:
            self.status_log.setText(f"[Ошибка 3D-рендера] Не удалось построить модель: {str(e)}")

    def on_search_error(self, error_message: str):
        self.status_log.setText(f"[База структур] {error_message}")

        self.set_placeholder_3d(error_message)
if __name__ == "__main__":
    sys.argv.append("--ignore-gpu-blocklist")
    sys.argv.append("--enable-webgl")
    sys.argv.append("--use-gl=desktop")
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
