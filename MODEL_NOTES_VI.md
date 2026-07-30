# Ghi chú mô hình và cách dùng

## 1. Luồng tính toán

`solve_schrodinger_poisson` tạo lưới z, xây dựng thế phi điều hòa, thế từ,
thế điện và thế Hartree. Ở mỗi vòng lặp:

1. giải bài toán trị riêng Hamiltonian sai phân hữu hạn;
2. chuẩn hóa hàm sóng theo `integral |psi|^2 dz = 1`;
3. giải mức Fermi từ điều kiện trung hòa điện tích;
4. tính mật độ electron;
5. giải Poisson;
6. trộn thế Hartree và kiểm tra đồng thời sai số thế lẫn sai số EF.

Sau khi hội tụ, `compute_transition_data` tính:

- chênh lệch năng lượng;
- phần tử ma trận lưỡng cực;
- form factor theo qz;
- tích phân bình phương form factor.

`compute_moap_direct` sau đó tính tổng theo q trong tọa độ trụ, tách:

- phonon quang / phonon áp điện;
- phát xạ / hấp thụ phonon;
- 1, 2, 3 photon.

## 2. Plot giống cách tổ chức bài [6]

- `plot_moap_spectrum`: tổng phổ và hai cơ chế phonon.
- `plot_moap_contributions`: từng bậc photon và từng quá trình phonon.
- `plot_sweep_spectra`: nhiều đường phổ khi thay B, T, Lz...
- `plot_linewidth_sweep`: FWHM theo tham số khảo sát.

## 3. Vì sao không dùng findpeaks

Các hàm tìm cực đại và FWHM được viết trực tiếp, nên không cần Signal
Processing Toolbox.

## 4. Hai kernel

- `direct_q_integral`: mặc định, ổn định hơn và tính trực tiếp trước phép
  khai triển Taylor-Maclaurin.
- `analytical_series`: đối chiếu PT. (5)/(7), dùng tổng hữu hạn theo s và
  eta; không nên dùng làm kết quả cuối nếu chưa kiểm tra hội tụ theo cả
  hai chỉ số.

## 5. Kiểm tra hội tụ cần thực hiện cho bài báo

Tăng lần lượt:

- `structure.Nz`;
- `oap.Nqz`, `oap.Nqperp`;
- `qz_max_inv_nm`, `qperp_max_inv_nm`;
- giảm bước năng lượng photon.

Kết quả được coi là hội tụ khi vị trí đỉnh, cường độ tương đối và FWHM
không còn thay đổi đáng kể.

## 6. Mật độ Debye mặc định

Mặc định dùng `n_e = 10^18 cm^-3` đúng theo bảng tham số trong file tính toán. Chế độ `from_sheet` vẫn được giữ để khảo sát, nhưng không phải mặc định.
