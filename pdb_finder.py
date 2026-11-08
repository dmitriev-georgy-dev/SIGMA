"""Модуль фонового поиска PDB-структур для SIGMA v1.1"""
import os
import hashlib
from PyQt6.QtCore import QObject, pyqtSignal
from protein_db import LOCAL_PROTEIN_DB

class LocalPDBFinder(QObject):
    """Класс для безопасного поиска файлов белков в локальной базе без использования интернета"""
    finished = pyqtSignal(str)  # Передаст путь к файлу при успехе
    error = pyqtSignal(str)     # Передаст текст ошибки, если белок не найден

    def __init__(self, protein_seq, resource_path_worker):
        super().__init__()
        self.protein_seq = protein_seq
        self.resource_path_worker = resource_path_worker

    def run(self):
        try:
            # Считаем MD5-хэш полученного после трансляции белка
            current_hash = hashlib.md5(self.protein_seq.encode()).hexdigest()

            if current_hash in LOCAL_PROTEIN_DB:
                file_name = LOCAL_PROTEIN_DB[current_hash]

                # Используем переданную функцию get_resource_path для совместимости с PyInstaller
                expected_path = self.resource_path_worker(os.path.join('structures', file_name))

                # Проверяем, физически ли лежит файл в папке
                if os.path.exists(expected_path):
                    self.finished.emit(expected_path)
                else:
                    self.error.emit(f"Файл {file_name} зарегистрирован в базе, но отсутствует в папке 'structures'.")
            else:
                self.error.emit("Данная аминокислотная последовательность отсутствует в локальной базе данных программы.")

        except Exception as e:
            self.error.emit(f"Внутренняя ошибка поиска: {str(e)}")
