# LUH3 to LUH1 Converter for TaiESM1/CLM4

## 概述

此工具將 LUH3 (CMIP7) 格式的土地利用資料轉換為 LUH1 格式，供 TaiESM1/CLM4 的 mksurfdat 工具使用。

## 轉換內容

### 輸入: LUH3 格式
- 解析度: 1440×720 (0.25°)
- PFT 結構: 15 natural PFTs + 64 CFTs (分離)
- 資料類型: double
- FillValue: -9999.0

### 輸出: LUH1 格式
- 解析度: 720×360 (0.5°)
- PFT 結構: 17 PFTs (混合)
- 資料類型: float
- FillValue: -999.0

## 主要變數對應

| LUH3 | LUH1 | 轉換邏輯 |
|------|------|----------|
| PCT_NAT_PFT[0-14] + PCT_NATVEG | PCT_PFT[0-14] | PCT_NATVEG × PCT_NAT_PFT[k] / 100 × residule / 100 |
| PCT_CROP | PCT_PFT[15] | residule - sum(PCT_PFT[0-14]) |
| N/A | PCT_PFT[16] | 0 (bare ground) |
| HARVEST_VH1/VH2/SH1/SH2/SH3 | HARVEST_VH1/VH2/SH1/SH2/SH3 | 空間降解析度 |
| GRAZING | GRAZING | 空間降解析度 |

## 安裝需求

```bash
# Python 3.6+
pip install numpy netCDF4 scipy
```

## 使用方法

### 1. 單一檔案轉換

```bash
python convert_luh3_to_luh1.py \
    -i /path/to/luh3_2019.nc \
    -o /path/to/luh1_2019.nc
```

### 2. 批次轉換 (1850-2023)

**步驟 1**: 編輯 `batch_convert_luh3_to_luh1.sh`，設定路徑：

```bash
LUH3_DIR="/your/luh3/data/directory"
LUH1_OUTPUT_DIR="/your/output/directory"
```

**步驟 2**: 調整檔案命名模式（如需要）：

```bash
# 例如，如果您的檔案名稱是:
# mksrf_landuse_clm6_histLUH3_1850.c251012.nc
# mksrf_landuse_clm6_histLUH3_1851.c251012.nc
# ...

LUH3_PATTERN="mksrf_landuse_clm6_histLUH3_YEAR.c251012.nc"
LUH1_PATTERN="mksrf_landuse_rcYEAR_converted.nc"
```

**步驟 3**: 執行批次轉換：

```bash
chmod +x batch_convert_luh3_to_luh1.sh
./batch_convert_luh3_to_luh1.sh
```

### 3. 驗證轉換結果

```bash
python validate_conversion.py \
    -i /path/to/luh1_2019.nc \
    --check-pft-sum \
    --check-ranges
```

## 轉換邏輯說明

### PFT 轉換

基於參考 R 腳本的邏輯：

```
1. 計算 residule = 100 - (glacier + lake + wetland + urban)

2. PFT 0-14 (自然植被):
   PCT_PFT[k] = PCT_NATVEG × PCT_NAT_PFT[k] / 100 × residule / 100

3. PFT 15 (作物):
   PCT_PFT[15] = residule - sum(PCT_PFT[0:15])
   
   邊界檢查:
   - 如果 < 0, 設為 0
   - 如果 sum > 100%, 調整為 100 - sum(PCT_PFT[0:15])

4. PFT 16 (裸地):
   PCT_PFT[16] = 0 (根據需求設定)

5. 重新縮放確保總和 = 100%:
   - 計算總和
   - 如果總和 ≠ 100%, 按比例縮放
   - 將剩餘誤差加到最大的 PFT
```

### 空間降解析度

使用 2×2 區塊平均法，從 0.25° 降到 0.5°：

```
對每個 LUH1 網格點 (i, j):
  取 LUH3 的 2×2 區塊 [2i:2i+2, 2j:2j+2]
  計算平均值 (忽略 NaN)
```

## 輸出檔案檢查

轉換完成後，請檢查：

1. **維度正確**
   ```bash
   ncdump -h output.nc | grep dimensions
   # 應該看到: lon = 720, lat = 360, pft = 17
   ```

2. **PFT 總和**
   ```bash
   # 使用 NCO 工具
   ncap2 -s 'pft_sum=PCT_PFT.total($pft)' output.nc test.nc
   ncwa -a pft -v PCT_PFT output.nc -O test.nc
   # 檢查 pft_sum 是否接近 100
   ```

3. **資料範圍合理**
   ```bash
   ncdump -v PCT_PFT output.nc | grep -A 5 "PCT_PFT ="
   # 檢查數值是否在 0-100 範圍內
   ```

## 已知問題與注意事項

### 1. Harvest 單位
- LUH3: `gC/m2/yr` (碳通量)
- LUH1: `unitless` (無單位)
- **當前版本**: 保持原始數值，未進行單位轉換
- **可能需要**: 根據 mksurfdat 的期望調整單位轉換係數

### 2. 2005-2006 交接處
- 如果您有 LUH1 原始資料 (1850-2005)，建議：
  - 1850-2005: 使用原始 LUH1
  - 2006-2023: 使用轉換後的 LUH3
  - 檢查 2005-2006 的連續性

### 3. 邊界像素
- 極地區域可能有較多 glacier/lake
- 請特別檢查高緯度地區的轉換結果

### 4. 記憶體需求
- 單一年份檔案: ~500 MB RAM
- 批次處理: 建議 2 GB+ RAM

## 後續步驟：使用 mksurfdat

轉換完成後，使用 TaiESM1 的 mksurfdat 工具：

```bash
# 範例 (實際指令請參考 TaiESM1 文件)
mksurfdat \
    --luh1-dir /path/to/converted/luh1/files \
    --start-year 1850 \
    --end-year 2023 \
    --output fsurdat.nc \
    --output-dyndat fdyndat.nc
```

## 問題排查

### 問題 1: "Failed to open NetCDF file"
- 檢查檔案路徑是否正確
- 確認有讀取權限
- 檢查檔案是否損壞: `ncdump -h file.nc`

### 問題 2: "PCT_PFT sum > 100%"
- 這是正常的中間狀態
- 腳本會自動重新縮放
- 如果最終輸出仍 > 100%, 請檢查 residule 計算

### 問題 3: 記憶體不足
- 減少批次處理的並行數
- 增加系統 swap 空間
- 使用分時段處理

## 參考資料

1. **R 腳本**: 原始 CLM4 PFT 映射邏輯
2. **C 程式**: CTSM5.2 land use data tool (百分比縮放邏輯)
3. **LUH3 文件**: https://luh.umd.edu/
4. **TaiESM1**: Taiwan Earth System Model

## 版本歷史

- v1.0 (2026-01-28): 初始版本
  - 基本 LUH3 → LUH1 轉換
  - PFT 結構轉換
  - 空間降解析度
  - Harvest 資料處理

## 授權

本工具基於 TaiESM1 計畫開發，用於科學研究目的。

## 聯絡資訊

如有問題或建議，請聯繫 TaiESM1 團隊。
# landuse_luh3_convert_luh1
