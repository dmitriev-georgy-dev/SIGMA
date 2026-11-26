"""Модуль фонового поиска PDB-структур для SIGMA v1.2"""
import os
import hashlib
from PyQt6.QtCore import QObject, pyqtSignal
from protein_db import PROTEIN_REGISTRY, LOCAL_PROTEIN_DB


class LocalPDBFinder(QObject):
    """
    Класс для безопасного поиска файлов белков в локальной базе без интернета.
    Поддерживает поиск по базовому шаблону (для точечных мутаций) и fallback по MD5-хэшу.
    """
    finished = pyqtSignal(str)  # Передаст путь к файлу при успехе
    error = pyqtSignal(str)     # Передаст текст ошибки, если белок не найден

    def __init__(self, protein_seq: str, resource_path_worker, base_protein_name: str = ""):
        super().__init__()
        self.protein_seq = protein_seq
        self.resource_path_worker = resource_path_worker
        self.base_protein_name = base_protein_name

    def run(self):
        try:
            file_name = None

            # 1. Приоритет: если выбран базовый белок из реестра (поддержка точечных мутаций)
            if self.base_protein_name and self.base_protein_name in PROTEIN_REGISTRY:
                file_name = PROTEIN_REGISTRY[self.base_protein_name]["file"]

            # 2. Запасной вариант: поиск точного совпадения по MD5-хэшу (для свободного ввода)
            if not file_name:
                current_hash = hashlib.md5(self.protein_seq.encode()).hexdigest()
                if current_hash in LOCAL_PROTEIN_DB:
                    file_name = LOCAL_PROTEIN_DB[current_hash]

            if file_name:
                expected_path = self.resource_path_worker(os.path.join('structures', file_name))

                if os.path.exists(expected_path):
                    self.finished.emit(expected_path)
                else:
                    self.error.emit(f"Файл {file_name} зарегистрирован в базе, но отсутствует в папке 'structures'.")
            else:
                self.error.emit("Данная последовательность отсутствует в локальной базе структур программы.")

        except Exception as e:
            self.error.emit(f"Внутренняя ошибка поиска: {str(e)}")
