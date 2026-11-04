"""Sequence Interface for Genetic Modeling and Analysis"""

import sys
import re
import os
from Bio.Seq import Seq
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QTextEdit, QVBoxLayout,
    QHBoxLayout, QWidget, QPushButton, QSplitter
)
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtCore import Qt, QThread, QObject, pyqtSignal

# Принудительно отключаем черные списки видеокарт Chromium на уровне ОС
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--ignore-gpu-blocklist --disable-gpu-driver-bug-workarounds"
os.environ["QT_ANGLE_PLATFORM"] = "d3d11"


def get_resource_path(relative_path):
    """Динамически определяет путь к ресурсам (работает и в .py, и в собранном .exe)"""
    if hasattr(sys, '_MEIPASS'):
        # Если запущено внутри собранного PyInstaller пакета
        return os.path.join(sys._MEIPASS, relative_path)
    # Если запущено как обычный скрипт .py
    return os.path.join(os.path.abspath("."), relative_path)


# Создаст папку 'structures' в директории скрипта, если её ещё нет (только при обычном запуске)
if not hasattr(sys, '_MEIPASS'):
    if not os.path.exists('structures'):
        os.makedirs('structures')

# Тестовые последовательности
Proinsuline = "tttgtgaaccaacacctgtgcggctcacacctggtggaagctctctacctagtgtgcggggaacgaggcttcttctacacacccaagacccgccgggaggcagaggacctgcaggtggggcaggtggagctgggcgggggccctggtgcaggcagcctgcagcccttggccctggaggggtccctgcagaagcgtggcattgtggaacaatgctgtaccagcatctgctccctctaccagctggagaactactgcaactag"
Lysozime = "aaagtctttggacgatgtgagctggcagcggctatgaagagtcacggacttgataactaccgcggatacagcctgggaaactgggtgtgtgctgcaaaattcgagagtaacttcaacacccaggctacaaaccgtaacaccgatggcagtaccgactacggcatcctacagatcaacagccgctggtggtgcaacgatggccgaaccccaggatcgaggaacctgtgcaacatcccgtgctcagccctgctgagctcagacataacagcgagcgtgaactgcgcgaagaagatcgtcagcgatggaaacggcatgaacgcgtgggtcgcctggcgcaaccgctgcaagggtaccgacgtccaggcgtggatcagaggctgccggctgtaa"
Myoglobin = "atgggtttgacggatggcgaatggcagctggtgctgaacgtgtggggcaaggtggaggccgacattccgggccacggccaggaggtgctgatccgcctgttcaagggccacccggagaccctggagaagttcgacaagttcaagcacctgaagagcgaggacgagatgaaggccagcgaggacctgaagaagcacggctgcaccgtgctgaccgccctgggcaccatcctgaagaagaagggccagcacgccgccgagatccagccgctggcccagagccacgccaccaagcacaagatcccggtgaagtacctggagttcatcagcgaggccatcatccaggtgctgcagagcaagcacccgggcgacttcggcgccgacgcccagggcgccatgaacaaggccctggagctgttccgcaaggacatcgccgccaagtacaaggagctgggctaccagggctag"


class LocalPDBFinder(QObject):
    """Класс для безопасного поиска файлов белков в локальной базе без использования интернета"""
    finished = pyqtSignal(str)  # Передаст путь к файлу при успехе
    error = pyqtSignal(str)     # Передаст текст ошибки, если белок не найден

    def __init__(self, protein_seq):
        super().__init__()
        self.protein_seq = protein_seq

    def run(self):
        try:
            import hashlib

            # Локальная база данных: хэш аминокислотной последовательности -> имя файла в папке structures
            LOCAL_PROTEIN_DB = {
                # Человеческий проинсулин (1zeg)
                hashlib.md5("FVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN".encode()).hexdigest(): "pdb1zeg.ent",
                # Лизоцим куриного яйца (1lyz)
                hashlib.md5("KVFGRCELAAAMKSHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL".encode()).hexdigest(): "pdb1lyz.ent",
                # Миоглобин (1mbo)
                hashlib.md5("MGLTDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGCTVLTALGTILKKKGQHAAEIQPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG".encode()).hexdigest(): "pdb1mbo.ent"
            }

            # Считаем MD5-хэш полученного после трансляции белка
            current_hash = hashlib.md5(self.protein_seq.encode()).hexdigest()

            if current_hash in LOCAL_PROTEIN_DB:
                file_name = LOCAL_PROTEIN_DB[current_hash]

                # КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: Обернули путь в функцию адаптации под PyInstaller
                expected_path = get_resource_path(os.path.join('structures', file_name))

                # Проверяем, физически ли лежит файл в папке
                if os.path.exists(expected_path):
                    self.finished.emit(expected_path)
                else:
                    self.error.emit(f"Файл {file_name} зарегистрирован в базе, но отсутствует в папке 'structures'.")
            else:
                self.error.emit("Данная аминокислотная последовательность отсутствует в локальной базе данных программы.")

        except Exception as e:
            self.error.emit(f"Внутренняя ошибка поиска: {str(e)}")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SIGMA")
        self.resize(1000, 600)
        self.setup_ui()

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Главный горизонтальный разделитель: Текст (слева) и Графика (справа)
        main_splitter = QSplitter(Qt.Orientation.Horizontal)

        # --- ЛЕВАЯ ЧАСТЬ (Работа с текстом) ---
        left_container = QWidget()
        left_layout = QVBoxLayout(left_container)
        left_layout.setContentsMargins(0, 0, 0, 0)

        # Вертикальный разделитель для редакторов ДНК и Белка
        text_splitter = QSplitter(Qt.Orientation.Vertical)

        self.dna_editor = QTextEdit()
        self.dna_editor.setPlaceholderText("Введите ДНК последовательность (A, T, G, C)...")
        self.dna_editor.setText(Proinsuline)

        self.protein_reader = QTextEdit()
        self.protein_reader.setPlaceholderText("Здесь появится белковая последовательность...")
        self.protein_reader.setReadOnly(True)

        text_splitter.addWidget(self.dna_editor)
        text_splitter.addWidget(self.protein_reader)

        self.protein_result_button = QPushButton("Транслировать в белок")
        self.protein_result_button.clicked.connect(self.translate_protein)

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
        main_splitter.setSizes([600, 400])

        final_layout = QHBoxLayout(main_widget)
        final_layout.addWidget(main_splitter)

    def on_search_success(self, file_path):
        self.protein_reader.append("\n[Система] Структура найдена в локальной базе.")
        self.render_protein_3d(file_path)

    def on_search_error(self, error_message):
        self.protein_reader.append(f"\n[Ошибка локальной базы] {error_message}")

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
            self.protein_reader.append("\n[Система] 3D-модель успешно отрисована!")
        except Exception as e:
            self.protein_reader.append(f"\n[Ошибка визуализации] Не удалось отобразить 3D: {str(e)}")

    def translate_protein(self):
        raw_dna = self.dna_editor.toPlainText().strip()
        dna_seq = re.sub(r'[\s\n\r]', '', raw_dna).upper()

        if not dna_seq:
            self.protein_reader.setText("Ошибка: Пожалуйста, введите последовательность ДНК.")
            return

        if not re.match(r'^[ATGC]+$', dna_seq):
            self.protein_reader.setText("Ошибка: Последовательность содержит недопустимые символы! Допускаются только A, T, G, C.")
            return

        if len(dna_seq) % 3 != 0:
            self.protein_reader.setText(f"Предупреждение: Длина последовательности ({len(dna_seq)}) не кратна 3. Возможен сдвиг рамки!")

        try:
            # Трасляция ДНК в белок с помощью Biopython
            rna = Seq(dna_seq).transcribe()
            protein = rna.translate(to_stop=True)
            self.protein_reader.setText(str(protein))

            # Создаем отдельный поток для локального поиска файла структуры
            self.thread = QThread()
            self.finder = LocalPDBFinder(str(protein))
            self.finder.moveToThread(self.thread)

            self.thread.started.connect(self.finder.run)
            self.finder.finished.connect(self.on_search_success)
            self.finder.error.connect(self.on_search_error)

            # Очистка памяти
            self.finder.finished.connect(self.thread.quit)
            self.finder.error.connect(self.thread.quit)
            self.finder.finished.connect(self.finder.deleteLater)
            self.finder.error.connect(self.finder.deleteLater)
            self.thread.finished.connect(self.thread.deleteLater)

            self.thread.start()
        except Exception as e:
            self.protein_reader.setText(f"Критическая ошибка обработки: {str(e)}")


if __name__ == "__main__":
    sys.argv.append("--ignore-gpu-blocklist")
    sys.argv.append("--enable-webgl")
    sys.argv.append("--use-gl=desktop")

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
