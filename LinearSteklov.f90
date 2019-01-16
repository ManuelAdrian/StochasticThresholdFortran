program LinearSteklov

    use mtmod
    use Parameters

    implicit none

    integer :: k

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T = 1000								!Tiempo de integración

    sigma_V  = 1.7_dp						!Intensidad de ruido asociado a los vectores
    sigma_H  = 1.7_dp						!Intensidad de ruido asociado a los humanos

	Lambda_V = 21000.0_dp					!Tasa de reclutamiento de los vectores
	beta_V 	 = 0.00003900042152404787_dp	!Tasa de contagio de los humanos hacia los vectores
	mu_V 	 = 2.1_dp						!Tasa de mortalidad de los vectores
	Nv		 = dble(Lambda_V)/dble(mu_V)

	Nh		 = 3650.0_dp
	beta_H 	 = 0.00003269533157348633_dp	!Tasa de contagio de los vectores hacia los humanos
    mu_H 	 = 0.0142857_dp					!Tasa de mortalidad de los humanos

	aReal8 	 = EPSILON(anyReal8)			!Valor real muy cercano a cero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Dt_reso = dble(T)/dble(N)	!Tamaño de paso de resolución
    Dt 		= R1*Dt_reso		!Tamaño de paso de la variable temporal

	x_zer0 = reshape((/ 2000.0_dp, 1.0_dp, 3500.0_dp, 150.0_dp /), (/ 1,4 /))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do k = 1,(L+1) 						!Vector del tiempo discretizado
		vect_time(k,1) = (k - 1)*Dt
	end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call Deterministic(x_zer0)
	call Stochastic_LS(x_zer0)

	call save_files(X_exact,X_stoch)
	
	contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine Deterministic(x_det)

	real(dp) :: x_det(1,4), F_matrix1(1,4)
	integer :: i

	X_exact(1,:) = x_det(1,:)

	do i = 1,L

		F_matrix1(1,1) = Lambda_V - beta_V*X_exact(i,1)*X_exact(i,4) - mu_V*X_exact(i,1)
		F_matrix1(1,2) = beta_V*X_exact(i,1)*X_exact(i,4) - mu_V*X_exact(i,2)
		F_matrix1(1,3) = mu_H*X_exact(i,4) - beta_H*X_exact(i,3)*X_exact(i,2)
		F_matrix1(1,4) = beta_H*X_exact(i,3)*X_exact(i,2) - mu_H*X_exact(i,4)

        X_exact(i+1,:) = X_exact(i,:) + F_matrix1(1,:)*Dt

	end do

	end subroutine Deterministic

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine Stochastic_LS(x_sto)

	real(dp) :: x_sto(1,4), Mat_a1(4,4), Mat_a2(4,4), G(2,4), dW(2,R1), Wdisc(1,2)
	real(dp) :: vec_a(4), vec_b(1,4), X_temp1(1,4), X_temp2(1,4), X_temp3(1,4)
	integer  :: seed, i, j, j1, j2
	logical  :: cond1, cond2, cond(4)

	seed = getseed()
	call sgrnd(seed)

	X_stoch(1,:) = x_sto(1,:)

	do i = 1,4					!CreaciÃ³n de las matrices a usar en el mÃ©todo LS
    	Mat_a1(i,:) = 0.0_dp
        Mat_a2(i,:) = 0.0_dp
        G(:,i) 	 	= 0.0_dp
    end do
 
    do j = 1,L
	    do j1 = 1,2
	        if (j == 1) then
	            dW(:,1) = 0.0_dp
                do j2 = 2,R1
	                dW(j1,j2) = sqrt(Dt_reso)*gaussrnd()
                end do
            else
                do j2 = 1,R1
	                dW(j1,j2) = sqrt(Dt_reso)*gaussrnd()
                end do
            end if
		end do

        vec_a(1) = - (beta_V*X_stoch(j,4) + mu_V)			!Numericamente mejor ponerlo así
        vec_a(2) = - mu_V
        vec_a(3) = - beta_H*X_stoch(j,2)
		vec_a(4) = - mu_H

        vec_b(1,1) = Lambda_V
        vec_b(1,2) = beta_V*X_stoch(j,4)*X_stoch(j,1)
        vec_b(1,3) = mu_H*X_stoch(j,4)
		vec_b(1,4) = beta_H*X_stoch(j,3)*X_stoch(j,2)

		do j1 = 1,4
			Mat_a1(j1,j1) = exp(Dt*vec_a(j1))		!Matriz A^1 en el mÃ©todo LS
			cond1 = DABS(vec_a(j1)) < aReal8
			cond2 = DABS(Mat_a1(j1,j1) - 1.0_dp) < aReal8

            if (cond1 .OR. cond2) then		!Matriz A^2 en el mÃ©todo LS.
				Mat_a2(j1,j1) = Dt
            else
				Mat_a2(j1,j1) = (dble(Mat_a1(j1,j1) - 1.0_dp)/dble(vec_a(j1)))
            end if
		end do

		G(1,1) = - sigma_V*X_stoch(j,4)*X_stoch(j,1)/dble(Nv)
        G(1,2) =   sigma_V*X_stoch(j,4)*X_stoch(j,1)/dble(Nv)
        G(2,3) = - sigma_H*X_stoch(j,2)*X_stoch(j,3)/dble(Nv)
        G(2,4) =   sigma_H*X_stoch(j,2)*X_stoch(j,3)/dble(Nv)
            
        Wdisc  = reshape((/ sum(dW(1,1:R1)), sum(dW(2,1:R1)) /),(/ 1,2 /))
            
        X_temp1(1,:) = matmul(X_stoch(j,:),Mat_a1)
        X_temp2 = matmul(vec_b,Mat_a2)
        X_temp3 = matmul(Wdisc,G)

        X_stoch(j+1,:) = X_temp1(1,:) + X_temp2(1,:) + X_temp3(1,:)

        do j1 = 1,4
			cond(j1) = X_stoch(j+1,j1) < 0
		end do

		if (cond(1) .OR. cond(2)) then
			X_stoch(j+1,1:2) = X_temp1(1,1:2) + X_temp2(1,1:2) - X_temp3(1,1:2)
		end if

		if (cond(3) .OR. cond(4)) then
			X_stoch(j+1,3:4) = X_temp1(1,3:4) + X_temp2(1,3:4) - X_temp3(1,3:4)
		end if

	end do

	end subroutine Stochastic_LS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine save_files(data_det,data_sto)

	real(dp) :: data_sto(L+1,4), data_det(L+1,4)
	integer  :: i, i1

100 format (1X,8000E17.6E3) 				!Printing format
    character(20) filename1, filename2

    filename1 = 'Mat_Det.dat'
    filename2 = 'Mat_Sto.dat'

    open(UNIT=45,FILE=filename1,position="APPEND")                                              
    open(UNIT=46,FILE=filename2,position="APPEND")

    do i = 1,(L+1)
		if (mod(i,100) == 0) then
        	write(45,100) vect_time(i,1), (data_det(i,i1),i1=1,4), (data_det(i,1) + data_det(i,2)),&
			& (data_det(i,3) + data_det(i,4))
        	write(46,100) vect_time(i,1), (data_sto(i,i1),i1=1,4), (data_sto(i,1) + data_sto(i,2)),&
			& (data_sto(i,3) + data_sto(i,4))
		end if
    end do
	
	close(45)
	close(46)

	end subroutine save_files

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end program LinearSteklov
