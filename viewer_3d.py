"""Модуль 3D-визуализации белков и подсветки мутаций для SIGMA v1.2"""
import warnings
import py3Dmol
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Protein3DViewer:
    """Генератор интерактивных 3D-сцен на базе py3Dmol с поддержкой маппинга остатков PDB"""

    def __init__(self):
        self.parser = PDBParser(QUIET=True)

    def _extract_residue_map(self, pdb_path: str) -> list[tuple[str, int]]:
        """
        Читает PDB-файл и возвращает упорядоченный список (chain_id, res_seq_num)
        для каждого аминокислотного остатка. Это связывает индекс в строке белка
        с реальным номером остатка внутри файла структуры.
        """
        residue_map = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = self.parser.get_structure("protein", pdb_path)

        # Берем первую модель (актуально для ЯМР-структур с несколькими конформациями)
        first_model = next(structure.get_models())
        for chain in first_model:
            for residue in chain:
                # hetflag == " " означает стандартный аминокислотный остаток (отсекаем воду HOH и лиганды)
                if residue.id[0] == " ":
                    residue_map.append((chain.id, residue.id[1]))

        return residue_map

    def build_protein_html(self, pdb_path: str, mutation_indices: list[int] = None) -> str:
        """
        Создает HTML-код с 3D-моделью белка.
        Если переданы mutation_indices, подсвечивает эти остатки красным и показывает боковые радикалы.
        """
        with open(pdb_path, "r", encoding="utf-8") as f:
            pdb_data = f.read()

        view = py3Dmol.view(width="100%", height="100%")
        view.addModel(pdb_data, "pdb")

        # Базовый стиль всей молекулы
        view.setStyle({"cartoon": {"color": "spectrum"}})

        # Подсветка точечных мутаций
        if mutation_indices:
            residue_map = self._extract_residue_map(pdb_path)

            for idx in mutation_indices:
                if 0 <= idx < len(residue_map):
                    chain_id, resi_num = residue_map[idx]
                    selection = {"chain": chain_id, "resi": resi_num}

                    # Выделяем мутировавший участок красным + рисуем боковой радикал (stick)
                    view.addStyle(
                        selection,
                        {
                            "cartoon": {"color": "#ff1744"},
                            "stick": {"colorscheme": "redCarbon", "radius": 0.22}
                        }
                    )

                    # Добавляем аккуратную текстовую метку позиции в 3D
                    view.addLabel(
                        f"Mut #{idx + 1}",
                        {
                            "position": selection,
                            "backgroundColor": "#b71c1c",
                            "fontColor": "white",
                            "fontSize": 11,
                            "backgroundOpacity": 0.85
                        },
                        selection
                    )

        view.zoomTo()
        return view._make_html()

    @staticmethod
    def build_blocked_html(reason: str) -> str:
        """Возвращает информативную заглушку, если 3D-отрисовка заблокирована анализатором"""
        return f"""
        <body style='background-color:#262626; color:#e0e0e0; display:flex; flex-direction:column;
                     justify-content:center; align-items:center; height:100vh; margin:0;
                     padding:30px; box-sizing:border-box; font-family:sans-serif; text-align:center;'>
            <div style='background-color:#332222; border:1px solid #d32f2f; border-radius:8px;
                        padding:24px; max-width:460px;'>
                <h3 style='color:#ff5252; margin-top:0;'>3D-визуализация недоступна</h3>
                <p style='font-size:14px; line-height:1.5; color:#cccccc; margin-bottom:0;'>
                    {reason}
                </p>
            </div>
        </body>
        """
