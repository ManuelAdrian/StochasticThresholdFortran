program Sde

    use mtmod
    use sde_par

    implicit none

    integer :: seed, i, i1, j, j1, j2, k, temp_var

100 format (1X,8000E17.6E3) 				!Printing format
    character(20) filename1, filename2, filename3, filename4, filename5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T = 10									!Tiempo de integración

    sigma_V  = 0.01_dp						!Intensidad de ruido asociado a los vectores
    sigma_H  = 0.01_dp						!Intensidad de ruido asociado a los humanos

    N_H 	 = 4000.0_dp					!Población de humanos constante

	Lambda_V = 60000.0_dp					!Tasa de reclutamiento de los vectores
	beta_V 	 = dble(0.438_dp)/dble(N_H)		!Tasa de contagio de los humanos hacia los vectores
	mu_V 	 = 2.1_dp						!Tasa de mortalidad de los vectores

	beta_H 	 = dble(0.01168_dp)/dble(N_H)	!Tasa de contagio de los vectores hacia los humanos
    mu_H 	 = 0.0142857_dp					!Tasa de mortalidad de los humanos

	aReal8 	 = EPSILON(anyReal8)			!Valor real muy cercano a cero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call cpu_time(start)

	seed = getseed()
	call sgrnd(seed)

    Dt_reso = dble(T)/dble(N)	!Tamaño de paso de resolución
    Dt 		= R1*Dt_reso		!Tamaño de paso de la variable temporal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nombrar y abrir archivos a utilizar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filename1 = 'Mat_Sv.dat'
    filename2 = 'Mat_Iv.dat'
    filename3 = 'Mat_Ih.dat'
    filename4 = 'extinction.dat'
    filename5 = 'time.dat'

    open(UNIT=45,FILE=filename1,position="APPEND")                                              
    open(UNIT=46,FILE=filename2,position="APPEND")
    open(UNIT=47,FILE=filename3,position="APPEND")
    open(UNIT=48,FILE=filename4,position="APPEND")

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vect_time(1,1) 		 = 0.0_dp	!Primera componente del vector de la variable temporal
    vect_time_print(1,1) = 0.0_dp	!Primera componente del vector a imprimir de la variable temporal
            
    do k = 1,L 						!Vector del tiempo discretizado
       vect_time(k + 1,1) = k*Dt
    end do

    do i = 1,N_simu					!Ciclo para el número total de simulaciones pedidas

10		k = 2
        temp_var = 0
        
        X = reshape((/ 1190.0_dp, 10.0_dp, 10.0_dp /), (/ 3,1 /))			!Condición inicial de SDE
	    X_exact = reshape((/ 1190.0_dp, 10.0_dp, 10.0_dp /), (/ 3,1 /))     !Condición inicial de ODE

        do i1 = 1,3					!Creación de las matrices a usar en el método LS
            Mat_a1(i1,:) = 0.0_dp
            Mat_a2(i1,:) = 0.0_dp
            G(i1,:) 	 = 0.0_dp
        end do
 
		if (i == 1) then			!Almacenar las condiciones iniciales de la ODE en una nueva matriz para guardarlas en un fichero
		    X1_exact(1,1) = X_exact(1,1)
    		X2_exact(1,1) = X_exact(2,1)
    		X3_exact(1,1) = X_exact(3,1)
		end if

!%%%%%%%% Almacenar las variables del SDEs en distintas matrices para guardarlas posteriormente en un fichero

       	Mat_Sv(1,i) = X(1,1)
       	Mat_Iv(1,i) = X(2,1)
       	Mat_Ih(1,i) = X(3,1) 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            vec_a(1) = - beta_V*X(3,1) - mu_V
            vec_a(2) = - mu_V
            vec_a(3) = - beta_H*X(2,1) - mu_H

            vec_b(1,1) = Lambda_V
            vec_b(2,1) = beta_V*X(1,1)*X(3,1)
            vec_b(3,1) = beta_H*N_H*X(2,1)

			do j1 = 1,3
				Mat_a1(j1,j1) = exp(Dt*vec_a(j1))		!Matriz A^1 en el método LS
                if (abs(vec_a(j1)) < aReal8) then		!Matriz A^2 en el método LS.
					Mat_a2(j1,j1) = Dt
                else
					Mat_a2(j1,j1) = (dble(Mat_a1(j1,j1) - 1.0_dp)/dble(vec_a(j1)))
                end if
			end do

            G(1,1) = - sigma_V*X(1,1)*X(3,1)
            G(2,1) =   sigma_V*X(1,1)*X(3,1)
            G(3,2) =   sigma_H*(N_H - X(3,1))*X(2,1)
            
            Wdisc  = reshape((/ sum(dW(1,1:R1)), sum(dW(2,1:R1)) /),(/ 2,1 /))

            X	= matmul(Mat_a1,X) + matmul(Mat_a2,vec_b) + matmul(G,Wdisc)

!			do j1 = 1,3
!				if (isnan(X(j1,1))) GOTO 10		!Si el programa bota basura, vuelve a mandar otra corrida
!			end do

            do j1 = 1,3
              if (X(j1,1) .GE. 0.0_dp) then		!Ya que se ha mostrado teoricamente la positividad de las soluciones, se procede a rectificar si hay un error numérico
                X(j1,1) =   X(j1,1)
                else
                X(j1,1) = - X(j1,1)
              end if
            end do

			if (i == 1) then					!Aproximación de la solución de la ODE
		    	vec_a1(1,1) = Lambda_V - beta_V*X_exact(1,1)*X_exact(3,1) - mu_V*X_exact(1,1)
		        vec_a1(2,1) = beta_V*X_exact(1,1)*X_exact(3,1) - mu_V*X_exact(2,1)
        		vec_a1(3,1) = beta_H*X_exact(2,1)*(N_H - X_exact(3,1)) - mu_H*X_exact(3,1)

		        X_exact = X_exact + vec_a1*Dt
            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			if (mod(j,T4) == 0) then			!Almacenar valores de las variables del SDE y ODE
            	vect_time_print(k,1) = vect_time(j+1,1)
            	Mat_Sv(k,i) = X(1,1) 
            	Mat_Iv(k,i) = X(2,1)
            	Mat_Ih(k,i) = X(3,1)

	            if (i == 1) then
			       	X1_exact(k,1) = X_exact(1,1)
       				X2_exact(k,1) = X_exact(2,1)
       				X3_exact(k,1) = X_exact(3,1)
                end if				
            	
				k = k + 1
 			end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
			if ((temp_var .LT. j) .AND. (X(2,1) .LT. 0.9999_dp) .AND.&
            	& (X(3,1) .LT. 0.9999_dp)) then			!Almacenar tiempo de extinción y valores de las variables de infección
               	
				time_ext = j + 1
   	            temp_var = L + 1

				temp_X(1) = X(2,1)						!Infectados vectores
                temp_X(2) = X(3,1)						!Infectados humanos
            end if
            
        end do

      	ext_time(i,1) = vect_time(time_ext,1)
      	ext_time(i,2) = temp_X(1)
      	ext_time(i,3) = temp_X(2)
    end do

    do i = 1,(L/T4)+1
        write(45,100) vect_time_print(i,1), X1_exact(i,1), (Mat_Sv(i,i1),i1=1,N_simu)
        write(46,100) vect_time_print(i,1), X2_exact(i,1), (Mat_Iv(i,i1),i1=1,N_simu)
        write(47,100) vect_time_print(i,1), X3_exact(i,1), (Mat_Ih(i,i1),i1=1,N_simu)
    end do

    do i = 1,N_simu
      write(48,100) (ext_time(i,j),j = 1,3)
    end do

    close(45)
    close(46)
    close(47)
    close(48)
    
	call cpu_time(finish)

	tim_fin = finish - start

    open(UNIT=49,FILE=filename5,position="APPEND")
    write(49,100) tim_fin, aReal8
    close(49)

end program Sde
