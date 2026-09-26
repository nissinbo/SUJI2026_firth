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

## CIでのPDF生成

GitHub Actionsの `.github/workflows/render-pdf.yml` は、保存済みの計算結果からスライドPDFだけを生成します。

- SASやRでモデルを再実行しません。
- `results/figures/` にコミット済みの図と `presentation/` の素材を読み込みます。
- Quarto 1.10.18を使い、Windows runner上で `scripts/render_slides.ps1` を実行します。
- PRと `main` へのpush、および手動実行で動きます。
- 生成した `presentation/slides.pdf` はActionsの `slides-pdf` artifactとして14日間保存します。
- `main` へのpush、または `main` を指定した手動実行では、生成・検証に成功したPDFを `github-actions[bot]` が `main` に自動コミットします。CI完了後に `git pull --ff-only` すると、手元のPDFも更新できます。
- PDFに差分がない場合やPRでは書き戻しません。生成中に `main` が先へ進んだ場合も古いPDFの書き戻しをスキップします。必要な場合は最新の `main` で手動実行してください。
- 書き戻し専用ジョブだけに `contents: write` を付与し、標準の `GITHUB_TOKEN` を使用します。この自動pushではCIを再起動しません。

したがって、計算結果を更新するときはローカル環境で結果・図を更新してコミットし、CIではrenderだけを再現します。
