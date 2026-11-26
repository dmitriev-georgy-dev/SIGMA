"""Модуль анализа генетических и белковых последовательностей для SIGMA v1.2"""
import re
from dataclasses import dataclass, field
from Bio.Seq import Seq


VALID_DNA_CHARS = set("ATGC")
VALID_PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class AnalysisResult:
    """Структурированный отчет о проверке и сравнении последовательностей"""
    is_valid: bool                          # Прошла ли последовательность базовую проверку символов
    error_message: str = ""                 # Текст ошибки, если ввод некорректен
    dna_seq: str = ""                       # Очищенная ДНК (если анализ шел от ДНК)
    protein_seq: str = ""                   # Итоговая белковая последовательность
    dna_mutations: list[int] = field(default_factory=list)      # Индексы измененных нуклеотидов
    protein_mutations: list[int] = field(default_factory=list)  # Индексы измененных аминокислот
    can_render_3d: bool = False             # Разрешено ли строить 3D-модель
    render_block_reason: str = ""           # Причина блокировки 3D (если can_render_3d == False)
    summary_log: str = ""                   # Готовая сводка для строки логов внизу окна


class SequenceAnalyzer:
    """Класс для валидации, трансляции и поиска точечных мутаций"""

    # Максимальная доля точечных замен (15%), при которой еще корректно опираться на исходный 3D-каркас
    MAX_MUTATION_RATIO = 0.15

    @staticmethod
    def clean_sequence(raw_text: str) -> str:
        """Убирает пробелы, переносы строк и переводит текст в верхний регистр"""
        return re.sub(r'[\s\n\r]', '', raw_text).upper()

    @staticmethod
    def find_invalid_chars(seq: str, allowed_chars: set[str]) -> list[str]:
        """Возвращает отсортированный список недопустимых символов в цепи"""
        return sorted(list(set(seq) - allowed_chars))

    @staticmethod
    def find_point_differences(ref_seq: str, cur_seq: str) -> list[int]:
        """Возвращает индексы позиций, в которых последовательности различаются (при равной длине)"""
        min_len = min(len(ref_seq), len(cur_seq))
        return [i for i in range(min_len) if ref_seq[i] != cur_seq[i]]

    def analyze_from_dna(self, raw_dna: str, ref_dna: str = "", ref_protein: str = "") -> AnalysisResult:
        """Полный цикл анализа при запуске от верхнего поля (ДНК -> Белок)"""
        dna_seq = self.clean_sequence(raw_dna)
        ref_dna_clean = self.clean_sequence(ref_dna)
        ref_prot_clean = self.clean_sequence(ref_protein)

        if not dna_seq:
            return AnalysisResult(is_valid=False, error_message="Поле ДНК пусто.")

        # 1. Проверка на недопустимые буквы в ДНК
        bad_chars = self.find_invalid_chars(dna_seq, VALID_DNA_CHARS)
        if bad_chars:
            return AnalysisResult(
                is_valid=False,
                error_message=f"Недопустимые символы в ДНК: {', '.join(bad_chars)}. Разрешены только A, T, G, C."
            )

        # 2. Трансляция ДНК в белок
        try:
            rna = Seq(dna_seq).transcribe()
            # Останавливаемся на первом стоп-кодоне
            protein_seq = str(rna.translate(to_stop=True))
        except Exception as e:
            return AnalysisResult(is_valid=False, error_message=f"Ошибка трансляции: {str(e)}")

        if not protein_seq:
            return AnalysisResult(
                is_valid=False,
                error_message="Синтез оборвался в самом начале (стоп-кодон в первом триплете)."
            )

        # 3. Если эталона нет (пользователь стер всё и ввел произвольную ДНК)
        if not ref_dna_clean or not ref_prot_clean:
            return AnalysisResult(
                is_valid=True,
                dna_seq=dna_seq,
                protein_seq=protein_seq,
                can_render_3d=True,
                summary_log="Трансляция произвольной ДНК завершена. Поиск точного совпадения в базе..."
            )

        # 4. Проверяем целостность рамки считывания и длины цепи
        if len(dna_seq) != len(ref_dna_clean):
            diff = len(dna_seq) - len(ref_dna_clean)
            action = f"вставлено {diff}" if diff > 0 else f"удалено {abs(diff)}"
            reason = f"Изменена длина ДНК ({action} нукл.). Сдвиг рамки считывания — 3D-моделирование заблокировано."
            return AnalysisResult(
                is_valid=True,
                dna_seq=dna_seq,
                protein_seq=protein_seq,
                can_render_3d=False,
                render_block_reason=reason,
                summary_log=reason
            )

        if len(protein_seq) < len(ref_prot_clean):
            reason = (
                f"Нонсенс-мутация: преждевременный стоп-кодон на позиции {len(protein_seq) + 1}. "
                "Цепь оборвана — 3D-моделирование заблокировано."
            )
            return AnalysisResult(
                is_valid=True,
                dna_seq=dna_seq,
                protein_seq=protein_seq,
                dna_mutations=self.find_point_differences(ref_dna_clean, dna_seq),
                can_render_3d=False,
                render_block_reason=reason,
                summary_log=reason
            )

        # 5. Поиск точечных мутаций (длина сохранена)
        dna_muts = self.find_point_differences(ref_dna_clean, dna_seq)
        prot_muts = self.find_point_differences(ref_prot_clean, protein_seq)

        # Проверяем, не заменено ли слишком много аминокислот
        if len(prot_muts) / len(ref_prot_clean) > self.MAX_MUTATION_RATIO:
            reason = "Слишком много замен в цепи (>15%). Исходный каркас недействителен — 3D заблокировано."
            return AnalysisResult(
                is_valid=True,
                dna_seq=dna_seq,
                protein_seq=protein_seq,
                dna_mutations=dna_muts,
                protein_mutations=prot_muts,
                can_render_3d=False,
                render_block_reason=reason,
                summary_log=reason
            )

        # Формируем понятный лог по результатам
        if not dna_muts:
            log = "Последовательность полностью совпадает с эталоном."
        elif not prot_muts:
            log = f"Обнаружено синонимичных замен в ДНК: {len(dna_muts)}. Структура белка не изменилась."
        else:
            log = f"Найдено точечных мутаций в ДНК: {len(dna_muts)}, замен аминокислот: {len(prot_muts)}."

        return AnalysisResult(
            is_valid=True,
            dna_seq=dna_seq,
            protein_seq=protein_seq,
            dna_mutations=dna_muts,
            protein_mutations=prot_muts,
            can_render_3d=True,
            summary_log=log
        )

    def analyze_from_protein(self, raw_protein: str, ref_protein: str = "") -> AnalysisResult:
        """Анализ при прямой правке белковой цепи в нижнем окне"""
        protein_seq = self.clean_sequence(raw_protein)
        ref_prot_clean = self.clean_sequence(ref_protein)

        if not protein_seq:
            return AnalysisResult(is_valid=False, error_message="Поле белковой последовательности пусто.")

        # 1. Проверка символов аминокислот (включая защиту от стоп-кодона '*' в середине)
        bad_chars = self.find_invalid_chars(protein_seq, VALID_PROTEIN_CHARS)
        if bad_chars:
            return AnalysisResult(
                is_valid=False,
                error_message=f"Недопустимые символы в белке: {', '.join(bad_chars)}. Используйте стандартные коды аминокислот."
            )

        if not ref_prot_clean:
            return AnalysisResult(
                is_valid=True,
                protein_seq=protein_seq,
                can_render_3d=True,
                summary_log="Проверка произвольной белковой цепи завершена."
            )

        # 2. Проверка длины белковой цепи
        if len(protein_seq) != len(ref_prot_clean):
            diff = len(protein_seq) - len(ref_prot_clean)
            action = f"добавлено {diff}" if diff > 0 else f"удалено {abs(diff)}"
            reason = f"Изменена длина белка ({action} аминокислот). Каркас нарушен — 3D-моделирование заблокировано."
            return AnalysisResult(
                is_valid=True,
                protein_seq=protein_seq,
                can_render_3d=False,
                render_block_reason=reason,
                summary_log=reason
            )

        # 3. Поиск точечных замен аминокислот
        prot_muts = self.find_point_differences(ref_prot_clean, protein_seq)

        if len(prot_muts) / len(ref_prot_clean) > self.MAX_MUTATION_RATIO:
            reason = "Слишком много точечных замен (>15%). 3D-моделирование заблокировано."
            return AnalysisResult(
                is_valid=True,
                protein_seq=protein_seq,
                protein_mutations=prot_muts,
                can_render_3d=False,
                render_block_reason=reason,
                summary_log=reason
            )

        log = (
            "Белковая цепь совпадает с эталоном."
            if not prot_muts
            else f"Обнаружено точечных замен аминокислот: {len(prot_muts)}."
        )

        return AnalysisResult(
            is_valid=True,
            protein_seq=protein_seq,
            protein_mutations=prot_muts,
            can_render_3d=True,
            summary_log=log
        )
