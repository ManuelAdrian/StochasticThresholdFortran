module Parameters

    implicit none

    integer, parameter :: dp = kind(0.d0)           ! double precision
    integer, parameter :: N = 1000000000, R1 = 100, L = N/R1

    integer :: T, Nh, Nv

    real(dp) :: Lambda_H, sigma_V, sigma_H, Lambda_V, beta_V, mu_V, beta_H, mu_H
    real(dp) :: Dt_reso, Dt, aReal8, anyReal8, start

    real(dp) :: x_zer0(1,4), vect_time(L+1,1), X_exact(L+1,4), X_stoch(L+1,4) 

end module Parameters
