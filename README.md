# SASとRにおけるFirth法の実装差の整理と数値的評価

SASユーザー総会2026（2026年10月15日）の発表資料と、再現用の解析コードを管理するリポジトリです。

- [発表スライド](presentation/slides.pdf)
- [発表原稿](presentation/speaker_notes_25min.md)
- [想定Q&A](presentation/qa.md)

## 構成

- `R/`: データ生成、モデル実行、集計、検証、作図
- `sas/`: PROC LOGISTIC・%FLによる解析
- `scripts/`: 実行・検証・スライド生成のエントリーポイント
- `data/`: 数値実験の入力とシナリオ定義
- `results/`: 保存済みのモデル出力、集計表、図、検証記録
- `presentation/`: Quartoスライド、発表原稿、Q&A

## 主な再現手順

プロジェクトルートから実行します。

```powershell
Rscript scripts/check_code.R
Rscript scripts/update_results.R
pwsh scripts/render_slides.ps1
```

`scripts/check_code.R` は保存済み結果との回帰チェックを行い、フルの数値実験は再実行しません。  
`scripts/update_results.R` は保存済みのモデル出力から集計表・図・検証記録を更新します。モデル当てはめや入力生成は別の明示的な処理です。

## 比較時の基本方針

同じFirth法でも、信頼区間・検定方式、標準誤差の定義、計算アルゴリズム、収束条件が実装ごとに異なります。本解析では、これらを分けて比較し、丸め前の数値と収束状態を確認します。
