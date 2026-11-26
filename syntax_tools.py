"""Модуль подсветки синтаксиса и точечных мутаций в текстовых полях для SIGMA v1.2"""
from PyQt6.QtGui import QSyntaxHighlighter, QTextCharFormat, QColor, QTextCursor, QFont
from PyQt6.QtWidgets import QTextEdit


class SequenceSyntaxHighlighter(QSyntaxHighlighter):
    """Подсветка недопустимых символов в реальном времени при вводе последовательности"""

    def __init__(self, document, allowed_chars: set[str]):
        super().__init__(document)
        self.allowed_chars = allowed_chars

        # Формат для ошибочных символов: красный фон + жирный белый шрифт
        self.error_format = QTextCharFormat()
        self.error_format.setBackground(QColor("#b71c1c"))
        self.error_format.setForeground(QColor("#ffffff"))
        self.error_format.setFontWeight(QFont.Weight.Bold)

    def highlightBlock(self, text: str):
        for idx, char in enumerate(text):
            # Игнорируем пробелы и переносы строк
            if char.isspace():
                continue
            if char.upper() not in self.allowed_chars:
                self.setFormat(idx, 1, self.error_format)


def clear_extra_highlights(text_edit: QTextEdit):
    """Сбрасывает всю ручную подсветку мутаций, оставляя стандартный цвет текста"""
    cursor = text_edit.textCursor()
    cursor.beginEditBlock()
    cursor.select(QTextCursor.SelectionType.Document)

    default_fmt = QTextCharFormat()
    cursor.setCharFormat(default_fmt)

    cursor.clearSelection()
    cursor.endEditBlock()


def highlight_mutations(text_edit: QTextEdit, mutation_indices: list[int], is_protein: bool = False):
    """
    Подсвечивает измененные символы по их индексам в очищенной последовательности.
    Корректно пропускает пробелы и переносы строк в самом QTextEdit.
    """
    clear_extra_highlights(text_edit)

    if not mutation_indices:
        return

    mut_set = set(mutation_indices)
    raw_text = text_edit.toPlainText()

    # Настраиваем цвет: для ДНК — янтарный, для белка — ярко-красный
    mut_format = QTextCharFormat()
    if is_protein:
        mut_format.setBackground(QColor("#d32f2f"))
        mut_format.setForeground(QColor("#ffffff"))
    else:
        mut_format.setBackground(QColor("#f57f17"))
        mut_format.setForeground(QColor("#000000"))
    mut_format.setFontWeight(QFont.Weight.Bold)

    cursor = text_edit.textCursor()
    cursor.beginEditBlock()

    clean_idx = 0
    for doc_pos, char in enumerate(raw_text):
        if char.isspace():
            continue

        if clean_idx in mut_set:
            cursor.setPosition(doc_pos)
            cursor.movePosition(QTextCursor.MoveOperation.Right, QTextCursor.MoveMode.KeepAnchor, 1)
            cursor.setCharFormat(mut_format)

        clean_idx += 1

    cursor.endEditBlock()
