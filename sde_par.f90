module sde_par

    implicit none

    integer, parameter :: dp = kind(0.d0)           ! double precision
    integer, parameter :: N_simu = 3, N = 10000000, R1 = 100, L = N/R1
    integer, parameter :: T4 = 10

    integer :: T, time_ext

    real(dp) :: N_H, sigma_V, sigma_H, Lambda_V, beta_V, mu_V, beta_H, mu_H
    real(dp) :: aReal8, anyReal8, start, finish, tim_fin
    real(dp) :: Dt_reso, Dt

    real(dp) :: Mat_a1(3,3), Mat_a2(3,3), X(3,1), vec_b(3,1)
    real(dp) :: vec_a(3), G(3,2), X_exact(3,1), vec_a1(3,1)

    real(dp) :: dW(2,R1), Wdisc(2,1), ext_time(N_simu,3)
    real(dp) :: vect_time(L+1,1), vect_time_print((L/T4)+1,1)

    real(dp) :: Mat_Ih((L/T4)+1,N_simu), Mat_Iv((L/T4)+1,N_simu), Mat_Sv((L/T4)+1,N_simu)
    real(dp) :: X1_exact((L/T4)+1,1), X2_exact((L/T4)+1,1), X3_exact((L/T4)+1,1), temp_X(2)

end module sde_par
