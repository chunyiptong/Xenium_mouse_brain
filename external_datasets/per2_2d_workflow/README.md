# Per2 smFISH：2D 操作流程

目標：每個偵測到的 spot 一列，包含原圖座標、原始強度量測、SCN ROI 和 cell ID。輸出是 Tab 分隔的 `.txt`，可直接交給 R。這份流程不執行 3D 分析。

另提供 [R 版本 `analyze_per2_2d.R`](analyze_per2_2d.R)，安裝與執行方式見 [README_R.md](README_R.md)。

## 1. 準備一張一致的 2D 影像

在 Fiji 以 Plugins → Bio-Formats → Bio-Formats Importer 開啟 ND2。

選擇其中一種方法，並對所有分析通道使用相同選擇：

- 單張切片：Image → Duplicate…，指定相同的一個 Z 層、需要的通道及時間點。檢查輸出只有一個 Z 層。
- 最大強度投影：Image → Stacks → Z Project… → Max Intensity，記錄開始與結束 Z 層。這是 2D 投影，會合併不同深度的訊號；強度和細胞歸屬也是投影上的量測。不要把這種表當作完整 3D 分子計數。

Image → Color → Split Channels，將 Per2 與 DAPI 各自用 File → Save As → Tiff… 儲存成單通道、單平面影像。勿轉成 RGB、8-bit、JPEG，也不要在分析圖上 Apply 亮度變更、畫入輪廓、裁切、縮放或旋轉。調整顯示對比供觀察可以，但不能改變原始強度。

在 Image → Properties… 記錄 Pixel width、Pixel height 和單位。腳本要求輸入 µm/pixel；若顯示為 pixel，請先從 ND2/顯微鏡 metadata 確認校正，勿自行猜測。

建立每張影像自己的資料夾，檔名例如：

```
sample01/
  Per2.tif
  DAPI.tif
  SCN_ROIs.zip
  DAPI_cp_masks.png
```

確認 Per2 和 DAPI 的通道身分，不能單憑顯示顏色推斷。原始 ND2 另外保留。

## 2. 畫 SCN contour

在相同像素尺寸與方向的 DAPI 或合併參考圖上，使用一致的解剖判定標準畫 SCN。若依靠標記通道，要避免依照 Per2 的強弱改變分析範圍。

1. Analyze → Tools → ROI Manager…。
2. 用 Polygon Selection 完成輪廓，按 `t`。
3. Rename… 為 `SCN_imageLeft`、`SCN_imageRight` 等唯一名稱。影像左右是否等於解剖左右要另外確認。
4. Deselect → More → Save…，存成 `SCN_ROIs.zip`。

腳本直接讀 ROI ZIP 的原始像素座標，不需要手動匯出 XY Coordinates。支援簡單 Polygon/Freehand/Traced ROI，不支援含洞或多區塊的 composite ROI。兩側 ROI 不應重疊；若重疊，腳本會標記 roi_overlap，並以分號列出 ROI 名稱。

## 3. 安裝分析工具（一次即可）

需要可用的 Python 3.10 以上版本及 Terminal。以下使用獨立環境，不修改系統 Python：

```bash
python3 -m venv "$HOME/per2-analysis-env"
source "$HOME/per2-analysis-env/bin/activate"
python -m pip install --upgrade pip
python -m pip install "cellpose[gui]" spotiflow numpy scipy tifffile pillow matplotlib roifile
```

第一次執行模型可能需要下載權重。之後每次開新 Terminal，先執行 `source "$HOME/per2-analysis-env/bin/activate"`。

## 4. 建立細胞或細胞核 label mask

```bash
python -m cellpose
```

將 DAPI.tif 拖入 Cellpose GUI，執行可用模型，檢查其結果是否真的是個別細胞核。不同版本的模型選項不同；不要因為模型能產生 mask 就假設它可靠。

逐一檢查並修正：黏在一起的核是否被拆開、同一個核是否被分成多塊、背景是否誤判。修正不佳的自動分割，必要時手繪。保留可編輯的 `_seg.npy`，並用 File 選單的「Save masks as PNG」（官方快捷鍵 Ctrl+N）另存整數 label mask，例如 `DAPI_cp_masks.png`。

正確 mask：背景值為 0；每個核／細胞各有不同整數值 1、2、3…。不要使用 screenshot、彩色 overlay、只有輪廓的影像，或所有細胞都等於 255 的二值圖。mask 與 Per2 的寬高必須完全一致。

細胞歸屬有三種明確模式：

- 只有 DAPI、只想可靠記錄核內 spots：`--mask-kind nuclei --expand-um 0`。輸出 assignment 為 `nucleus_only`，表中的 cell_id 此時只是核所代表的細胞識別碼，不表示測到全細胞。
- 只有 DAPI、接受估計核周區域：`--mask-kind nuclei --expand-um 距離`。腳本按物理距離有限擴張核標籤，重疊區按最近核邊界分配；這些是估計的領域，不是真實細胞膜。等距邊界由演算法決定，應檢查邊界標記。距離必須以代表影像和解剖知識驗證，沒有通用值。
- 有足夠的全細胞染色、已檢查完整細胞 mask：`--mask-kind cells --expand-um 0`。

只有 DAPI 時，不可能保證所有胞質 spots 的真實細胞歸屬。無歸屬的點會保留，不能強迫每個 spot 都有 cell ID。2D 投影中重疊的細胞也可能無法區分。

## 5. 執行 Spotiflow、強度量測和 cell assignment

先只處理一張測試影像。以下 `cd` 是本工作區的實際路徑；`sample01` 需自己建立並放入上面四個檔案。

```bash
cd /Users/chunyiptong/Downloads/20261010/per2_2d_workflow
source "$HOME/per2-analysis-env/bin/activate"
python analyze_per2_2d.py \
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
  --out sample01/results_run01
```

必須把 `YOUR_PIXEL_SIZE_X`、`YOUR_PIXEL_SIZE_Y` 改成實測的 µm/pixel 數字；把 source-description 改成實際資料。若使用最大投影，改為 `--image-mode max_projection`。預設命令只指派核內 spots；若採用經驗證的核周擴張距離，把 `--expand-um 0` 改成該 µm 值。

腳本以 Spotiflow `general` 模型偵測 2D spots；預設使用模型的偵測閾值。若需改變，可加 `--prob-thresh 0.數值`（範圍 0–1）。先人工驗證代表影像，再鎖定同一分析策略。調整參數後使用新的輸出資料夾；腳本拒絕覆寫非空目錄。

## 6. 輸出與量測定義

- `per2_spots.txt`：所有偵測點，包含 SCN 內外及未指派點。
- `cells.txt`：每個 label 一列，包含零個 spots 的物件。提供全部點數與 SCN 內點數。roi_at_centroid 是「質心所在的 ROI」，不是整個細胞都在 SCN 內的保證。
- `contours.txt`：有順序的 ROI 頂點，單位 pixel。
- `assignment_labels.tif`：實際用來指派的完整 label mask。
- `qc_overlay.png`：Per2 背景、SCN 輪廓、cell/nucleus 邊界及 spots 疊圖。
- `run_info.json`：輸入參數、版本、座標及強度定義。原始 ROI、mask、ND2、程式也要保留。

位置：`x_px, y_px` 是原圖 0-based 像素座標，x 向右、y 向下。`x_um, y_um` 是相對原圖原點的 µm 座標，不是顯微鏡 stage 座標。小數座標保留，mask lookup 使用最近像素。

強度：由原始輸入的 2D Per2 圖量測，沒有用模型分數代替強度。

- peak_raw：以 spot 座標為中心、半徑 2 pixels 的圓內最高強度。
- background_local：半徑 3–5 pixels 圓環的強度中位數，排除其他已偵測 spots 的量測圓。少於 10 個可用背景像素時輸出 NaN。
- integrated_bgsub：量測圓內的強度總和，減去 background_local × 像素數。可能為負，不強制歸零。
- score：Spotiflow 的偵測分數，與 intensity 不同。

上述 2/3/5 pixels 只是起始參數，必須根據解析度及 spot 尺寸驗證，可透過 `--spot-radius`、`--bg-inner`、`--bg-outer` 修改。不同解析度的影像需考慮一致物理尺寸的量測區域。量測圓未限制在 cell mask 內，避免 cell 邊界截斷量測；但因此附近其他細胞的訊號也可能影響強度。

此版本沒有做重疊 spots 的訊號解混。aperture_overlaps_other_spot、near_mask_boundary、edge_window、roi_overlap 是人工 QC 用的標記，不是已自動解決的問題。顯微鏡的實際飽和上限也要另外檢查（儲存成 16-bit 不代表訊號真的有 16-bit 動態範圍）。

## 7. 人工檢查後才批次

先挑至少幾張強、弱及背景較高的代表影像，檢查：

1. SCN contour 和 mask 是否與 Per2 完全對齊。
2. Spotiflow 是否漏掉弱點、誤判背景、將亮團錯誤計數。
3. 單細胞／單核分割是否正確。
4. 核周估計是否侵入鄰居；未指派比例是否合理。
5. 原圖量測區域是否包含完整 spot，背景圓環是否受到密集訊號污染。

overview PNG 只供快速檢查；密集區需放大原圖，必要時將 assignment_labels.tif 和 spots 座標載入 napari 檢查。核分割或參數不理想時，修正並重跑。完成後再用同一規則處理其他影像。

## 8. 在 R 讀取

```r
result_dir <- "/Users/chunyiptong/Downloads/20261010/per2_2d_workflow/sample01/results_run01"
spots <- read.delim(file.path(result_dir, "per2_spots.txt"),
                    na.strings = c("", "NA", "nan"), stringsAsFactors = FALSE)
cells <- read.delim(file.path(result_dir, "cells.txt"),
                    na.strings = c("", "NA", "nan"), stringsAsFactors = FALSE)

# 所有 SCN 內 spots，包括尚未指派細胞的點
scn_spots <- spots[!is.na(spots$roi_id), ]

# 已指派的 SCN spots；查看 assignment 欄，分辨 nuclei/cell/estimated。
assigned <- scn_spots[!is.na(scn_spots$cell_id), ]

# cells 保留零個 spots 的細胞／核，可直接使用 n_spots_in_scn。
head(assigned)
head(cells)
```

正式統計前合併 animal_id、CT、section_id、染色與拍攝批次等 metadata。若只測核內 spots，不應將 n_spots_in_scn 解讀為全細胞轉錄本數。偵測點也未必一對一等同單分子，尤其重疊 spots、亮團或轉錄位點。

## 驗證範圍

程式已用合成資料檢查 ROI 配對、座標尺度、已知強度／背景、mask 指派、有限距離擴張、零 spots 細胞、空偵測表及完整檔案輸出。尚未執行真實 Per2 影像的 Spotiflow 推論，亦未驗證你的影像上的 Cellpose 品質、模型閾值或量測參數。

## 官方參考

- Fiji 投影：https://imagej.net/imaging/z-functions
- Cellpose GUI：https://cellpose.readthedocs.io/en/latest/gui.html
- Cellpose mask 格式：https://cellpose.readthedocs.io/en/latest/outputs.html
- Spotiflow 安裝與 2D API：https://github.com/weigertlab/spotiflow
- Spotiflow API：https://weigertlab.org/spotiflow/api.html
