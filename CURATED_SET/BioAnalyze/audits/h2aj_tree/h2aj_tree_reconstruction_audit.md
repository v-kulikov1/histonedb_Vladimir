# Аудит реконструкции пайплайна дерева H2A.J / cH2A / nonplacental

## Что восстановлено

Восстановлен воспроизводимый CLI:

- `CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py`

Он делает три вещи:

1. Копирует минимально нужные historical-артефакты в external storage:
   `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\tree`
2. Воспроизводит `historical` и `clean` наборы входов для дерева в:
   `CURATED_SET/BioAnalyze/data/h2aj_tree`
3. Собирает evidence-папку и воспроизводит historical rooting для нуклеотидного дерева.

Основной результат этой реконструкции:

- подтверждённый успешный web-run в PhyML найден только для белкового дерева от `3 июня 2025`;
- mammal/nonplacental-белковая итерация восстановлена по сохранённым FASTA/PHY, но её web-result bundle не найден;
- провал нуклеотидного дерева объясняется прежде всего upstream-подготовкой последовательностей:
  BLAST-фрагменты вместо полного CDS, обрезание хвоста, фильтрация `SQK/SHK` по выровненной позиции, смешение белковых и нуклеотидных входов, а затем ручная чистка коротких записей.

## Приоритет источников

При конфликте источников использовался такой порядок доверия:

1. архивные outputs и bundle-файлы:
   `All_AA_0206.phy`, `All_AA_1006.phy`, `All_NUC_1606.phy`, `All_NUC_0907.phy`, `all_aa_0206_1__phy_phyml.zip`, `nuc_tree_bootstrap.nwk`, `rooted_tree.nwk`
2. сохранённые FASTA-файлы из `Гистоны` и `My_test`
3. Jupyter notebooks из `My_test` и `Работа над грантом HistoneJ`
4. Telegram export
5. `Draft2.docx`

Это важно, потому что память из переписки в одном месте говорит про `HKY85`, а реальный PhyML bundle показывает другой запуск:

- `aa`
- модель `LG`
- `BioNJ`
- `FreeRate(3)`
- `SH-like branch supports`
- без bootstrap

Следовательно, bundle сильнее любой позднейшей verbal-memory.

## Ключевая хронология

### 15 мая 2025

В notebook `15_05 BLAST N .ipynb` нуклеотидные кандидаты собираются как BLAST-HSP-фрагменты:

- используется `hsp.sbjct.replace('-', '')`

То есть берётся не гарантированно полный CDS, а subject-кусок из BLAST-выравнивания.

В той же линии видно и раннюю трансляцию проблемных нуклеотидов:

- есть предупреждение про `partial codon`

Это первый жёсткий маркер того, что pipeline уже на входе работал не с полными CDS.

### До 2 июня 2025

`Draft2.docx`, особенно слайды `4-7`, подтверждает early-концепцию работы:

- сначала BLAST по человеческому H2A.J;
- дальше берутся записи с `"J"` в title;
- затем делается alignment;
- затем последовательности делятся на `SQK` и `SHK` по позиции `132`;
- потом строится дерево;
- после этого интерпретируется, что mammalian `SQK` образует отдельную ветвь.

Это очень важный источник: он показывает, что исходная логика поиска была motif/title-driven, а не CDS-grounded.

### 2 июня 2025

Сохранены следующие файлы:

- `H2AJ.fasta` = `230` записей
- `cH2A.fasta` = `57` записей
- `All_AA_0206.fasta` = `387` записей
- `All_AA_0206.phy` = `387 208`

Критическая поправка к памяти:

- `230 + 57 + 51 = 338`, а не `387`

Значит `All_AA_0206` не равен позднему mammal/nonplacental-набору. Это ранний broad-outgroup белковый запуск, что хорошо согласуется с notebook `02_06 File_for_alignment_outgroup.ipynb`.

### 3 июня 2025

`all_aa_0206_1__phy_phyml.zip` подтверждает реальный успешный запуск на сайте PhyML:

- дата web-exec в stdout соответствует `3 июня 2025`
- data type = `aa`
- model = `LG`
- starting tree = `BioNJ`
- rate model = `FreeRate`
- classes = `3`
- supports = `SH-like branch supports`
- taxa = `387`

Это и есть strongest confirmed tree-build event.

### 10-11 июня 2025

Следующая белковая итерация уже соответствует той памяти, где используются:

- H2A.J
- mammalian cH2A
- nonplacental H2A

Сохранилось:

- `All_AA_1006.fasta`
- `All_AA_1006.phy`

Именно здесь сходится арифметика:

- `H2AJ = 230`
- `cH2A = 57`
- `platypus/nonplacental = 51`
- сумма = `338`

`All_AA_1006.phy` начинается с:

- `338 180`

То есть mammal/nonplacental-белковый набор был, но напрямую подтверждённого PhyML bundle для него не найдено.

### 16 июня 2025

Сохранился первый большой нуклеотидный набор:

- `Nuc_merged_sequences.fasta`
- `All_NUC_1606.fasta`
- `All_NUC_1606.phy`

Контрольная точка:

- `All_NUC_1606.phy = 315 396`

Notebook `Важные скрипты.ipynb` подтверждает, как для этой стадии вытягивались `cH2A_nuc` через `coded_by` из protein-accession.

### 17 июня 2025

Notebook `02_06 File_for_alignment_outgroup.ipynb` показывает точный historical rooting workflow:

1. читается `nuc_tree_bootstrap.nwk`
2. ищутся все листья с `cH2A_nuc`
3. проверяется `is_monophyletic(...)`
4. проверка падает
5. затем используется fallback с одной `cH2A_nuc`

Сохранённый `rooted_tree.nwk` совпадает с автоматически восстановленным вариантом, где outgroup = первая найденная `cH2A_nuc`.

### 7 июля 2025

В Telegram есть критическая самодиагностика:

- в “нуклеотидное” дерево попали белковые H2A.J вместо нуклеотидных

Это уже не косвенный, а прямой источник причины сбоя.

### 8 июля 2025

Сохраняются:

- `SQK_nuc.fasta` = `227`
- `protein_from_SQK_nuc.fasta` = `227`
- `SQK_nuc(without short).fasta` = `224`
- `protein_from_SQK_nuc(without short).fasta` = `224`

То есть к этому моменту pipeline уже работает в режиме manual rescue:

- сначала берутся нуклеотиды,
- потом они переводятся в белок,
- потом вручную отбрасывается всё, что “после перевода не выглядит как H2A.J / не содержит SQK”.

Это очень сильный индикатор того, что исходный нуклеотидный набор был нестабилен.

### 9 июля 2025

Появляется merged-набор:

- `cH2A_nuc + platypus_nuc + SQK_nuc.fasta` = `312`
- `All_NUC_0907.phy` = `312 396`

Это уже июльская пересборка после ручной чистки.

### 18 июля 2025

В переписке фиксируется влияние short sequences:

- короткие H2A.J-последовательности ложатся в сторону outgroup/canonical area

Это согласуется с `without short`-веткой от `8 июля 2025`.

### 23 июля 2025

В переписке отмечается, что нуклеотидное дерево всё ещё “странное”:

- часть marsupials падает в H2A.J-clade
- часть H2A.J уезжает отдельно

То есть даже после ручной чистки upstream-проблема не исчезла.

### 30 июля 2025

Это финальная и самая важная диагностика причины:

- использовались BLAST `aligned sequences`
- из-за этого последовательности обрезались
- C-tail с `SQK/SHK` часто не попадал целиком

Именно это лучше всего связывает вместе BLAST-HSP pipeline, motif-driven filtering и провал нуклеотидного дерева.

## Восстановленный ход работы

### 1. Ранняя логика поиска H2A.J

Самый ранний слой работы был не “сначала полный CDS, потом ортологизация”, а скорее:

- взять человеческий H2A.J
- сделать BLAST
- оставить то, что похоже по title/annotation
- выровнять
- посмотреть мотив в хвосте
- разделить на `SQK` и `SHK`
- интерпретировать дерево

Это видно и по `Draft2.docx`, и по `Копия Фильтрация SHK_SQK.ipynb`.

### 2. Белковая ветка

Белковая ветка распадается на две historical-итерации.

Итерация A:

- `All_AA_0206 = 387 / 208`
- early broad-outgroup protein tree
- именно она подтверждена архивным PhyML bundle

Итерация B:

- `All_AA_1006 = 338 / 180`
- H2A.J + mammalian cH2A + nonplacental set
- исторически nonplacental-набор назывался `platypus`

### 3. Нуклеотидная ветка

Нуклеотидная ветка шла в несколько этапов:

1. BLAST-derived fragments
2. очистка и перевод в белок
3. motif-based фильтрация `SQK/SHK`
4. сборка `All_NUC_1606`
5. rooting attempt с полным `cH2A_nuc`
6. ручная июльская пересборка `SQK_nuc`
7. удаление short sequences
8. merge в `312` записей

Это не один чистый pipeline, а серия исправлений после того, как выяснялось, что предыдущий набор собран не так.

## Почему не получилось построить нормальное нуклеотидное дерево

### Причина 1. Искались не полные CDS, а BLAST-HSP фрагменты

Самая ранняя проблема зафиксирована в `15_05 BLAST N .ipynb`:

- используется `hsp.sbjct.replace('-', '')`

Это значит:

- в FASTA уходил ровно кусок hit-а из BLAST alignment,
- а не обязательный полный coding sequence.

Если цель различать `SQK` и `SHK` по C-tail, такой pipeline заранее опасен, потому что BLAST-hit может не дотягиваться до конца CDS.

### Причина 2. Эти фрагменты потом переводились в белок

В той же линии работы есть предупреждение:

- `partial codon`

То есть часть последовательностей уже была не кратна трём. Если после этого их переводить и потом принимать решения по белковому мотиву, классификация становится шумной.

### Причина 3. Фильтрация `SQK/SHK` была привязана к выровненной позиции, а не к надёжно восстановленному хвосту

`Копия Фильтрация SHK_SQK.ipynb` показывает позиционную логику:

- `target_position = 131`

То есть motif проверялся на фиксированной aligned-позиции. Но если исходная последовательность уже обрезана по BLAST fragment, такая aligned-position логика может “видеть” не настоящий C-tail, а смещённый или неполный остаток.

### Причина 4. В нуклеотидное дерево попали белковые H2A.J

В переписке от `7 июля 2025` прямо сказано:

- были вставлены белковые последовательности H2A.J вместо нуклеотидных

Это уже само по себе достаточно, чтобы сломать interpretation tree-input.

### Причина 5. Затем пришлось вручную спасать набор через перевод и motif-cleanup

Вместо одного стабильного входа появляются:

- `SQK_nuc`
- `protein_from_SQK_nuc`
- `SQK_nuc(without short)`
- `protein_from_SQK_nuc(without short)`

То есть pipeline перешёл в режим:

- взять нуклеотиды,
- перевести в белок,
- проверить, это всё ещё H2A.J или нет,
- выкинуть короткие / подозрительные записи.

Такой rescue-stage полезен как диагностика, но он подтверждает, что исходная нуклеотидная выборка была собрана нестабильно.

### Причина 6. Short sequences дополнительно искажали topology

Отдельный июльский вывод:

- короткие записи тянули часть H2A.J к canonical/outgroup области

Это объясняет, почему даже после ручной фильтрации topology оставалась “странной”.

### Причина 7. Даже rooting-этап был нестабилен, потому что весь `cH2A_nuc` не монофилетичен

Это не первичная причина wrong input, но важное дополнительное осложнение.

Historical notebook и автоматически восстановленный `postprocess` показали:

- весь набор `cH2A_nuc` нельзя безопасно взять как outgroup;
- MRCA этого множества уже включает огромный кусок H2A.J-последовательностей;
- поэтому использовался fallback на одной `cH2A_nuc`.

Иными словами, даже после сборки дерева рутование по “всем canonical H2A_nuc” не работало.

## Главный вывод

Нуклеотидное дерево не “сломалось из-за плохой модели PhyML”.

Оно ломалось раньше:

1. при поиске кандидатов через BLAST aligned fragments
2. при обрезании хвоста
3. при переводе неполных фрагментов
4. при motif-фильтрации по fixed aligned position
5. при смешении белковых и нуклеотидных входов

Самая компактная формулировка причины:

- исходно SQK-кандидаты искались и обрезались так, что C-конец, на котором держится различение `SQK/SHK`, часто был неполным или смещённым

После этого downstream-дерево уже не могло быть устойчивым.

## Важные поправки к historical memory

### Поправка 1

`387` не равно `H2AJ + cH2A + nonplacental`.

Правильная интерпретация:

- `387` относится к early broad-outgroup белковому запуску `All_AA_0206`
- `338 = 230 + 57 + 51` относится к later mammal/nonplacental белковому набору `All_AA_1006`

### Поправка 2

Подтверждённый успешный web-run, найденный в архиве, был не `HKY85`, а:

- `LG`
- `BioNJ`
- `FreeRate(3)`
- `SH-like branch supports`

### Поправка 3

Historical label `platypus` в этих файлах фактически обозначает broader nonplacental-набор, а не только утконоса.

Поэтому clean-выходы нормализуют это в `nonplacental`.

## Что лежит в восстановленных выходах

### External archive

- `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\tree\tree_archive_manifest.tsv`

### Repo outputs

- `CURATED_SET/BioAnalyze/data/h2aj_tree/evidence/`
- `CURATED_SET/BioAnalyze/data/h2aj_tree/historical/`
- `CURATED_SET/BioAnalyze/data/h2aj_tree/clean/`

Особенно полезные файлы:

- `evidence/phyml_aa_0206_runbook.md`
- `evidence/notebook_evidence.tsv`
- `evidence/chat_evidence.tsv`
- `clean/postprocess/nuc_tree_all_cH2A_outgroup_failure.txt`

## Итог

Скрипт и аудит восстанавливают не только “какие файлы были поданы в PhyML”, но и сам historical failure mode:

- белковое дерево получилось на early broad-outgroup наборе;
- later mammal/nonplacental белковый набор был подготовлен отдельно;
- нуклеотидная ветка поехала из-за того, как искались SQK-кандидаты и как потом обрезались BLAST-derived последовательности;
- это дополнительно ухудшилось смешением модальностей и короткими последовательностями;
- rooting по полному `cH2A_nuc` тоже не работал, потому что outgroup был немонофилетичен.
## 2026-04-01 clean override

The `clean` outputs intentionally diverge from `historical` in one extra way:

- they drop three short H2A.J/SQK records:
- `Myotis-lucifugus|XM_006084274`
- `Homo-sapiens|AK303301`
- `Homo-sapiens|AL133626`
- July clean SQK outputs now use `SQK_nuc(without short)` and `protein_from_SQK_nuc(without short)` as the canonical H2A.J source.
- `clean/postprocess/*.nwk` remain historical reference trees and must not be treated as trees inferred from the newly filtered clean alignments.
