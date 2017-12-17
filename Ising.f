      program ISING_2D_metropolis
        !Simulation of 2D Ising model with zero external field using Metropolis algorithm
        !Bibliography:
        !**************************************************************
        !Toral, R., & Colet, P. (2014). Stochastic numerical methods:
        !an introduction for students and scientists. John Wiley & Sons.
        !Chapter 5
        !
        !Tomás Sintes' notes:
        !https://ifisc.uib-csic.es/users/tomas/MFCS/SquareLatticeIsing.pdf
        !**************************************************************

        !REMARKS:
        !-->Position (x,y) in grid is labbeled with with integer i=(y-1)*L+x
        !-->To compute the effect of correlations between measures in the error
        ! we have used an aproximation (R. Toral notes)
        implicit none
        integer*4 L,N
        parameter (L=50, N=L*L) !Lattice dimensions
        double precision T,mag,sus,Cv,E !Observables
        double precision tau,cm,cE,err_m,err_E !To compute errors
        double precision hamil,dran_u !Real functions
        double precision h(-4:4)
        real start,finish !To compute CPU time
        integer*4 i_dran,sum_n !Integer functions
        integer*4 neigh(4,N),s(N)
        double precision aux,last_value,last_value_E !Real dummy variables
        integer*4 i,j,k,m,steps,step_therm !Integer dummy variables
        call dran_ini(1994)
        open(unit=1, file="mag_vs_T.dat", status="unknown")
        open(unit=2, file="E_vs_T.dat", status="unknown")
        neigh=0
        call Create_neigbours(neigh,L,N)
        call random_IC(s,N) !random initial condition
        call cpu_time(start)
        do k=1,10
        T=3.0d0-0.15d0*dble(k) !temperature in J/Kb units

        !Compute possible acceptance h(i)=min(1,exp(-beta*DH({s})))
        !As there are only 4 NN and spin values are -1 or +1--> 5 different values of h
        !h(-4),h(-2),h(0),h(2),h(4)
          do i=-4,0,2
            h(i)=1
          enddo
          h(2)=dexp(-4.0d0/T)
          h(4)=dexp(-8.0d0/T)

          !THERMALIZATION
          steps=1000*N !1000 MC steps
          do i=1,steps
            call rejection(N,s,neigh,h)
          enddo

          !INITIALIZE MEASURES
          mag=0.0d0 !Magnetization
          sus=0.0d0 !Magnetic susceptibility
          E=0.0d0 !Energy
          Cv=0.0d0 !Specific heat
          cm=0.0d0 !one step correlation function (magnetization)
          cE=0.0d0
          tau=0.0d0
          last_value=abs(dble(sum(s)))
          last_value_E=hamil(s,N,T,neigh)
          !MAIN PROGRAM
          !+++++++++++++++++++++++++++++++++++++++++++++++++++
          steps=1000 !Number of measures
          do i=1,steps
            ! This loop pretends to decrease the correlation between measures
            do j=1,N  !ONE mc step
              call rejection(N,s,neigh,h) !PROPOSAL/REJECTION
            enddo

            !MEASURE
            aux=abs(dble(sum(s)))
            mag=mag+aux
            sus=sus+aux**2.0d0
            cm=cm+aux*last_value
            last_value=aux

            aux=hamil(s,N,T,neigh)
            E=E+aux !It would be much more efficient to update E at each step
            Cv=Cv+aux**2.0d0
            cE=cE+aux*last_value_E
            last_value_E=aux

          enddo
          mag=mag/(dble(steps)*dble(N)) !Magnetization is in [0,1]
          sus=(sus/(dble(steps)*dble(N)**2.0d0)-mag**2.0d0)
          cm=(cm/(dble(steps)*dble(N)**2.0d0)-mag**2.0d0)/sus
          if (cm.ne.1.0) then !if cm=1 <m^2>=<m(i)m(i+1)>
            tau=cm/(1.0d0-cm)
          endif
          write(*,*)"----------------------------------------"
          write(*,*)sus,cm,tau
          err_m=sqrt(sus*(2*tau+1)/dble(steps))
          sus=sus/T

          E=E/dble(steps)
          Cv=(Cv/dble(steps)-E**2.0d0)
          cE=(cE/dble(steps)-E**2.0d0)/Cv
          if (cE.ne.1.0) then !if cE=1 <E^2>=<E(i)E(i+1)>
            tau=cE/(1.0d0-cE)
          endif
          write(*,*)cV,cE,tau
          err_E=sqrt(cV*(2*tau+1)/dble(steps))
          Cv=Cv/T**2.0d0

          write(1,*) T,mag,err_m
          write(2,*) T,E,err_E

        enddo
        call cpu_time(finish)
        write(*,*) 'Time in seconds',finish-start
        close(1)
        close(2)
      end

      subroutine rejection(N,s,neigh,h)
        integer*4 N,rd,prod_spin,sum_n,i_dran
        integer*4 s(N),neigh(4,N)
        double precision dran_u,u,h(-4:4)
        !MAKES REJECTION STEP
        rd=i_dran(N) !choose spin at random
        prod_spin=s(rd)*sum_n(rd,N,s,neigh)
        if (prod_spin.le.0) then
        !This first if avoid to invoke a pseudo random number at each step
        !if prod_spin<=0 h=1--> change is always accepted
          s(rd)=-s(rd)
        else if (h(prod_spin).ge.dran_u()) then
          s(rd)=-s(rd)
        endif
        return
      end

      integer*4 function sum_n(rd,N,s,neigh)
        integer*4 rd,N,i
        integer*4 s(N),neigh(4,N)
        !sum spins of all nearest neighbours of site rd
        sum_n=0
        do i=1,4
          sum_n=sum_n+s(neigh(i,rd))
        enddo
        return
      end

      double precision function hamil(s,N,T,neigh)
        !compute Hamiltonian.
        !Actually, I'm computing Hamiltonian/J
        double precision T
        integer*4 N
        integer*4 s(N),neigh(4,N)
        integer*4 i,j
        hamil=0.0d0
        do i=1,N
          do j=1,4 !4=number of NN
            hamil=hamil-s(i)*s(neigh(j,i))
          enddo
        enddo
        return
      end function

      subroutine random_IC(s,N)
        integer*4 i,N,rd
        integer*4 s(N)
        !RANDOM INITIAL CONDITION
        do i=1,N
          rd=i_dran(2)
          s(i)=2*rd-3 !random number, 1 or -1
        enddo
        return
      end

      subroutine Create_neigbours(neigh,L,N)
        !neigh(i,j)= neighbour i of site j
        ! i=1--> right neighbour
        ! i=2--> left neighbour
        ! i=3--> down neighbour
        ! i=4--> up neighbour
        !The implementations of the (periodic) Boundary Conditions is included in this array
        integer*4 k,m,N,L
        integer*4 neigh(4,N)


        do k=1,L !coordinate x (coulmns)
          do m=1,L !coordinate y (rows)
            j=(m-1)*L+k

            ! right neighbours
            neigh(1,j)=j+1
            if (k.eq.L) then
              neigh(1,j)=(m-1)*L+1 !first coulmn
            endif

            !left neighbours
            neigh(2,j)=j-1
            if (k.eq.1) then
              neigh(2,j)=(m-1)*L+L !last coulmn
            endif

            !down neighbours
            neigh(3,j)=j-L
            if (m.eq.1) then
              neigh(3,j)=k+L*(L-1) !last row
            endif

            !up neighbpurs
            neigh(4,j)=j+L
            if (m.eq.L) then
              neigh(4,j)=k !first row BC
            endif

          enddo
        enddo

        return
      end subroutine
