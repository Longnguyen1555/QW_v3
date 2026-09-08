# Hướng dẫn chạy nhanh

## Bước 1: kiểm tra folder

Mở MATLAB, đặt **Current Folder** tại thư mục `QW_MOAP_DeltaDoped`, sau đó chạy:

```matlab
Main_quick_test
```

Kết quả đúng phải kết thúc bằng dòng:

```text
[QUICK TEST] PASS
```

## Bước 2: chạy mô phỏng đầy đủ

```matlab
Main
```

Các file `.mat`, `.csv`, `.png` được ghi vào `results/`.

## Bước 3: chọn loại kết quả

Trong `Main.m`, sửa:

```matlab
cfg.run.task = 'single';
```

Các lựa chọn:

- `single`: một cấu hình và một phổ;
- `sweep_spectra`: nhiều phổ khi quét một tham số;
- `linewidth_sweep`: FWHM/HWHM theo tham số;
- `all_demo`: chạy cả ba nhóm kết quả.

## Bước 4: quét B, T hoặc Lz

Ví dụ quét từ trường:

```matlab
cfg.sweep.parameter = 'B_T';
cfg.sweep.values = [5 10 15];
```

Quét nhiệt độ:

```matlab
cfg.sweep.parameter = 'T_K';
cfg.sweep.values = [77 150 300];
```

Quét độ rộng đặc trưng của giếng:

```matlab
cfg.sweep.parameter = 'Lz_nm';
cfg.sweep.values = [4.5 5.0 5.5];
```

## Bước 5: đổi cơ chế phonon hoặc bậc photon

```matlab
cfg.oap.mechanisms = {'optical', 'piezoelectric'};
cfg.oap.photon_orders = [1 2];
```

Có thể chọn riêng một cơ chế:

```matlab
cfg.oap.mechanisms = {'optical'};
```

## Bước 6: chạy toàn bộ kiểm thử MATLAB

```matlab
addpath(genpath(pwd));
run_all_tests
```
