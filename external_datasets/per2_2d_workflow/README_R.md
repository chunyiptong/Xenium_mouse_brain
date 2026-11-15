# Per2 2D：R 版本

`analyze_per2_2d.R` 是 `analyze_per2_2d.py` 的 R 版本。保留原有命令列參數、六個輸出檔名及表格欄位。影像讀取、Fiji ROI 解析、強度量測、細胞／細胞核歸屬、核周擴張及 QC 圖由 R 執行。

Spotiflow 偵測沿用 Python 的 `general` 模型，透過 `reticulate` 呼叫。若使用 `--points-tsv` 提供已偵測的座標和分數，整個流程不需要 Python。程式不執行 Cellpose 分割、通道選取或 Z 投影；輸入準備與量測定義見 [README.md](README.md)。

## 在 Terminal 偵測，再交給 R

### 原生 RData 輸出

`spotiflow_to_rdata.py` 將結果存成 R 可直接 `load()` 的 `results.RData`，並另存 `points.tsv`。在執行 Spotiflow 的 Python 環境中安裝依賴：

```bash
python -m pip install spotiflow numpy tifffile rdata

python spotiflow_to_rdata.py \
  --img "~/Downloads/20261010/vip-cck-per2-cor/CT2/49B_Per2.tif"
```

省略 `--out` 時，會在原圖旁建立 `49B_Per2_spotiflow_rdata` 資料夾；也可用 `--out "/完整路徑/新結果資料夾"` 指定。已有非空資料夾時拒絕覆寫。

`--prob-thresh 0.2` 可設定探索用的偵測閾值；省略時使用模型預設閾值，並非取消篩選。0.2 僅為設定方式示範，需依原圖與人工標記驗證。

```r
result_dir <- path.expand("~/Downloads/20261010/vip-cck-per2-cor/CT2/49B_Per2_spotiflow_rdata")
load(file.path(result_dir, "results.RData"))

# points 是原始 N × 2 座標矩陣，兩欄為 y、x，均為 0-based pixels。
colnames(points) <- c("y_px", "x_px")
head(points)

# details 是具名 list，保留所有模型回傳欄位及陣列軸序。
details$prob
details$intens
dim(details$heatmap)
dim(details$flow)

# 若需要含座標、分數和強度的表格：
spots <- read.delim(file.path(result_dir, "points.tsv"), na.strings = "nan")
```

讀取 `.RData` 不需要 R 額外套件或 Python。NumPy 浮點陣列轉為 R 的 double；不改變陣列維度，也不增加原始資料的精度。`None` 轉為 `NULL`，`NaN`／`Inf` 保留。`load()` 會建立或取代目前環境內同名的 `points`、`details`。

`tests/check_spotiflow_rdata.py` 已以合成資料驗證 R 原生讀取、座標與陣列軸序、空結果、巢狀欄位、NaN／Inf／NULL、TSV、輸入檢查與拒絕覆寫；尚未執行真實 Spotiflow 模型推論。

### JSON 輸出

`spotiflow_to_r.py` 可獨立執行 Spotiflow，將完整結果匯出給 R，不需要 `reticulate`，也不要求 R 和 Python 使用相同 CPU 架構：

```bash
python spotiflow_to_r.py --img "/完整路徑/Per2.tif" --out "/完整路徑/spotiflow_results"
```

在已安裝 `spotiflow`、`numpy`、`tifffile` 的 Python 環境執行。輸入須為單平面灰階 2D TIFF。省略 `--out` 時，輸出到原圖旁的 `<影像檔名>_spotiflow` 資料夾；已有非空資料夾時拒絕覆寫。

- `points.tsv`：`y_px`、`x_px`、`score`（`details.prob`）、`intensity`（`details.intens`）。座標仍為 0-based。這個檔案也可用於 `analyze_per2_2d.R --points-tsv`。
- `details.json.gz`：保留 `details` 的所有欄位，包含完整 heatmap、flow、subpix 與逐點分數／強度。陣列維度順序不變；無法提供的欄位與非有限數值以 JSON `null` 儲存。壓縮減少檔案大小，R 讀取後仍需容納完整陣列的記憶體。

```r
source("/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/read_spotiflow.R")
result <- read_spotiflow("/完整路徑/spotiflow_results")
points <- result$points
details <- result$details
head(points)
dim(details$heatmap)
dim(details$flow)
```

`intensity` 沿用 Spotiflow 回傳的 `details.intens`，不是原分析流程中扣除背景的 `integrated_bgsub`。格式依據 [Spotiflow 官方 API](https://weigertlab.org/spotiflow/api.html)。匯出／R 讀取已用合成輸出驗證，包括陣列軸序、空偵測及缺失值；尚未執行真實模型推論。

## 安裝

在 R／RStudio 執行一次：

```r
install.packages(c("ijtiff", "png", "jsonlite"))

# 只有需要由 R 啟動 Spotiflow 偵測時才需要：
install.packages("reticulate")
```

使用 TIFF／PNG 即可；若讀取灰階 JPEG，另外安裝 `jpeg`。定量影像仍應使用原始、無損的 TIFF，避免 JPEG 壓縮改變強度。

需要 Spotiflow 時，先依原流程安裝 Python 環境與 `spotiflow`，並在啟動 R 前設定：

```bash
export RETICULATE_PYTHON="$HOME/per2-analysis-env/bin/python"
```

也可在尚未載入 Python 的 R 工作階段中設定：

```r
Sys.setenv(RETICULATE_PYTHON = path.expand("~/per2-analysis-env/bin/python"))
```

R 與 Python 必須使用相同 CPU 架構。本機檢查時 R 為 `x86_64`，Homebrew Python 為 `arm64`；使用 `reticulate` 推論前需使用架構一致的 R／Python 環境。已有座標表的模式不受這項限制。腳本不會自動安裝 Python 套件或更換環境。

## 從 Terminal 執行

原來的命令只需將 `python analyze_per2_2d.py` 改為 `Rscript analyze_per2_2d.R`，其餘參數相同。以下路徑與校正值必須替換為自己的資料：

```bash
cd /Users/chunyiptong/Downloads/20261010/per2_2d_workflow

Rscript analyze_per2_2d.R \
  --per2 sample01/Per2.tif \
  --mask sample01/DAPI_cp_masks.png \
  --mask-kind nuclei \
  --expand-um 0 \
  --rois sample01/SCN_ROIs.zip \
  --image-id sample01 \
  --image-mode single_plane \
  --source-description "填入原始ND2檔名、Per2通道、Z層或投影範圍、時間點" \
  --pixel-size-x YOUR_PIXEL_SIZE_X \
  --pixel-size-y YOUR_PIXEL_SIZE_Y \
  --out sample01/results_R_run01
```

`YOUR_PIXEL_SIZE_X`、`YOUR_PIXEL_SIZE_Y` 是 µm/pixel；程式不會猜測校正值。其他選項：

| 參數 | 預設／用途 |
| --- | --- |
| `--mask-kind` | 必填：`nuclei` 或 `cells` |
| `--image-mode` | 必填：`single_plane` 或 `max_projection`，僅記錄輸入類型 |
| `--expand-um` | `0`；僅供 `nuclei` 使用，依 X/Y 實際像素尺寸擴張 |
| `--spot-radius` | `2` pixels |
| `--bg-inner` | `3` pixels |
| `--bg-outer` | `5` pixels |
| `--prob-thresh` | 省略時使用 Spotiflow 模型預設值 |
| `--points-tsv` | 已有的座標／分數表；提供時跳過 Spotiflow |
| `--help` | 顯示用法 |

輸出目錄必須不存在或為空；程式會拒絕覆寫非空目錄。

## 已有 spots 座標：純 R 分析

`points.tsv` 為 Tab 分隔，至少有三欄：

```text
y_px	x_px	score
120.25	310.50	0.95
132.00	325.75	0.88
```

座標必須是原圖的 **0-based 像素座標**，不是 R 的 1-based 矩陣索引，也不是 µm。可直接使用原 Python 輸出的 `per2_spots.txt` 作為 `--points-tsv`；其他欄位會被忽略。只有表頭、沒有資料列也可以。

在上述命令加入 `--points-tsv sample01/points.tsv` 即可。

在 RStudio 可先 `source()` 再呼叫 `main()`；`source()` 本身不會開始分析：

```r
source("/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/analyze_per2_2d.R")

pixel_size_x <- NA_real_  # 改成實測的 µm/pixel
pixel_size_y <- NA_real_  # 改成實測的 µm/pixel
sample_dir <- "/完整路徑/sample01"

result <- main(c(
  "--per2", file.path(sample_dir, "Per2.tif"),
  "--mask", file.path(sample_dir, "DAPI_cp_masks.png"),
  "--mask-kind", "nuclei",
  "--rois", file.path(sample_dir, "SCN_ROIs.zip"),
  "--image-id", "sample01",
  "--image-mode", "single_plane",
  "--source-description", "填入原始ND2檔名、通道、Z層、時間點",
  "--pixel-size-x", as.character(pixel_size_x),
  "--pixel-size-y", as.character(pixel_size_y),
  "--points-tsv", file.path(sample_dir, "points.tsv"),
  "--out", file.path(sample_dir, "results_R_run01")
))

head(result$spots)
head(result$cells)
```

## 輸出與相容性

輸出為 `per2_spots.txt`、`cells.txt`、`contours.txt`、`assignment_labels.tif`、`qc_overlay.png`、`run_info.json`。前三項是 Tab 分隔表，R 版將字串與表頭加上引號，讀取後欄位內容與 Python 版相同。布林值使用 `True`／`False`，缺失的量測值使用 `nan`。原 README 的 `read.delim()` 範例可直接使用。

- 原始 TIFF 強度不縮放；PNG 讀入後還原原始 bit depth，保留 16-bit 強度及 label ID。
- 輸出 label mask 使用 32-bit unsigned integer TIFF。
- X/Y、ROI 頂點及質心均使用 0-based；查詢 mask 時才轉為 R 索引。ROI 包含邊界點，支援整數與小數頂點、ROI ZIP 和中文內嵌名稱。
- 核周擴張採 Euclidean distance，使用各自的 X/Y µm/pixel；等距時按來源像素的 X、Y 升序決定。
- 強度、背景及重疊處理沿用原版；不去除負值，也不做重疊訊號解混。
- QC 圖沿用原版顏色與方向，字型、座標軸與版面由 R 繪圖系統呈現。
- `run_info.json` 額外記錄 R 與 R 套件版本，使用模型時記錄 Python 模型套件版本。
- 讀入範圍為灰階 TIFF／PNG（原圖另支援灰階 JPEG）。彩色圖、多通道／多平面圖與 composite ROI 會被拒絕。R 版另要求輸入分數為有限數值。

## 驗證

`tests/compare_r_python.py` 使用相同合成資料分別執行兩個版本，逐欄比較三張表、逐像素比較 assignment mask，並驗證 32-bit TIFF、QC PNG、metadata 及拒絕覆寫。涵蓋：

- 16-bit TIFF／PNG、32-bit label、大於 255 的 label ID。
- 細胞、細胞核、非等向性像素與有限距離擴張。
- 整數／小數 spots、影像邊緣、重疊 spots、重疊 ROI、ROI 邊界。
- 小數 ROI 頂點、中文 ROI 名稱、單一 ROI 與 ROI ZIP。
- 零 spots 的細胞、空偵測表、已知原圖強度及背景。
- 24 組隨機距離轉換，包含等距選擇、單列與單行影像。

執行驗證需要原 Python 程式的依賴及上述 R 套件，但不需要 Spotiflow：

```bash
python tests/compare_r_python.py
```

尚未在本機執行真實 Spotiflow 模型推論；測試不代表已驗證實際影像上的偵測品質或量測參數。

套件與格式參考：[ijtiff 讀取](https://docs.ropensci.org/ijtiff/reference/read_tif.html)、[ijtiff 寫入](https://docs.ropensci.org/ijtiff/reference/write_tif.html)、[ImageJ ROI 格式](https://github.com/imagej/ImageJ/blob/master/ij/io/RoiDecoder.java)。
