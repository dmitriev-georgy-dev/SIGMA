"""Sequence Interface for Genetic Modeling and Analysis (SIGMA) v1.1"""

import sys
import re
import os
from Bio.Seq import Seq
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QTextEdit, QVBoxLayout,
    QHBoxLayout, QWidget, QPushButton, QSplitter, QComboBox, QLabel
)
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtCore import Qt, QThread

# Импортируем поисковый движок и реестр базы данных
from pdb_finder import LocalPDBFinder
from protein_db import PROTEIN_REGISTRY

# Принудительно отключаем черные списки видеокарт Chromium на уровне ОС
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--ignore-gpu-blocklist --disable-gpu-driver-bug-workarounds"
os.environ["QT_ANGLE_PLATFORM"] = "d3d11"


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
        self.resize(1100, 650)
        self.setup_ui()

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Главный вертикальный макет, чтобы в самый низ положить строку логов
        window_layout = QVBoxLayout(main_widget)
        window_layout.setContentsMargins(5, 5, 0, 5)

        # Основной горизонтальный разделитель (Текст / Графика)
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # --- ЛЕВАЯ ЧАСТЬ (Работа с текстом) ---
        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        left_layout.setContentsMargins(0, 0, 0, 0)

        # Новый элемент: Выпадающий список выбора белков
        self.protein_selector = QComboBox()
        self.protein_selector.addItems(PROTEIN_REGISTRY.keys())
        self.protein_selector.currentIndexChanged.connect(self.on_protein_selected)

        text_splitter = QSplitter(Qt.Orientation.Vertical)

        self.dna_editor = QTextEdit()
        self.dna_editor.setPlaceholderText("Здесь отобразится кодирующая ДНК выбранного белка...")
        self.dna_editor.setReadOnly(True)  # Теперь пользователь не редактирует вручную

        self.protein_reader = QTextEdit()
        self.protein_reader.setPlaceholderText("Здесь появится белковая последовательность...")
        self.protein_reader.setReadOnly(True)

        text_splitter.addWidget(self.dna_editor)
        text_splitter.addWidget(self.protein_reader)

        self.protein_result_button = QPushButton("Транслировать и смоделировать 3D")
        self.protein_result_button.clicked.connect(self.translate_protein)

        left_layout.addWidget(QLabel("Выберите белок из базы данных:"))
        left_layout.addWidget(self.protein_selector)
        left_layout.addWidget(text_splitter)
        left_layout.addWidget(self.protein_result_button)

        # --- ПРАВАЯ ЧАСТЬ (Интерактивная 3D Графика) ---
        self.graphics_panel = QWebEngineView()
        self.graphics_panel.setHtml(
            "<body style='background-color:#2e2e2e; color:#888; display:flex; justify-content:center; "
            "align-items:center; height:100vh; font-family:sans-serif;'>Ожидание структуры белка...</body>"
        )

        main_splitter.addWidget(left_container)
        main_splitter.addWidget(self.graphics_panel)
        main_splitter.setSizes([550, 550])

        # --- НИЖНЯЯ ЧАСТЬ (Строка логов) ---
        self.status_log = QLabel("[Система] Программа готова к работе. Выберите биополимер.")
        self.status_log.setStyleSheet("color: #aaa; background-color: #252525; padding: 5px; font-family: monospace;")

        # Собираем всё в главное окно
        window_layout.addWidget(main_splitter)
        window_layout.addWidget(self.status_log)

        # Инициализируем первое значение
        self.on_protein_selected()

    def on_protein_selected(self):
        """Автоматически подставляет ДНК при смене выбранного белка в списке"""
        selected_name = self.protein_selector.currentText()
        if selected_name in PROTEIN_REGISTRY:
            self.dna_editor.setText(PROTEIN_REGISTRY[selected_name]["dna"])
            self.protein_reader.clear()
            self.status_log.setText(f"[Система] Загружена последовательность ДНК для: {selected_name}")

    def on_search_success(self, file_path):
        self.status_log.setText("[Система] Структура успешно найдена в локальной базе. Начинается отрисовка...")
        self.render_protein_3d(file_path)

    def on_search_error(self, error_message):
        self.status_log.setText(f"[Ошибка базы] {error_message}")

    def render_protein_3d(self, file_path):
        import py3Dmol
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                pdb_data = f.read()

            view = py3Dmol.view(width="100%", height="100%")
            view.addModel(pdb_data, "pdb")
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.zoomTo()

            html_content = view._make_html()
            self.graphics_panel.setHtml(html_content)
            self.status_log.setText("[Система] 3D-модель биополимера успешно отрисована!")
        except Exception as e:
            self.status_log.setText(f"[Ошибка визуализации] Не удалось отобразить 3D: {str(e)}")

    def translate_protein(self):
        raw_dna = self.dna_editor.toPlainText().strip()
        dna_seq = re.sub(r'[\s\n\r]', '', raw_dna).upper()

        if not dna_seq:
            self.status_log.setText("[Предупреждение] Поле ДНК пусто.")
            return

        try:
            # Трансляция ДНК в белок
            rna = Seq(dna_seq).transcribe()
            protein = rna.translate(to_stop=True)
            self.protein_reader.setText(str(protein))
            self.status_log.setText("[Система] Трансляция завершена. Поиск структуры в потоке...")

            # Поиск структуры через QThread
            self.thread = QThread()
            self.finder = LocalPDBFinder(str(protein), get_resource_path)
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
        except Exception as e:
            self.status_log.setText(f"[Критическая ошибка] Ошибка обработки: {str(e)}")


if __name__ == "__main__":
    sys.argv.append("--ignore-gpu-blocklist")
    sys.argv.append("--enable-webgl")
    sys.argv.append("--use-gl=desktop")

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
