
module simulation_utils
    use omp_lib
    implicit none
contains

    subroutine initialize_neighbors(LEFT_WALL,RIGHT_WALL,L,NEIGHBORS)

    implicit none
    integer :: ij
    integer, intent(in) :: LEFT_WALL(:), RIGHT_WALL(:)
    integer,intent(inout) :: NEIGHBORS(:,:)
    integer, intent(in) :: L

      do ij = 1, L**2

        if (any(RIGHT_WALL == ij)) then
            NEIGHBORS(ij,1) = ij - L +1
        else
            NEIGHBORS(ij,1) = ij +1
        end if

        if (any(LEFT_WALL == ij)) then
            NEIGHBORS(ij,2) = ij +L -1
        else
            NEIGHBORS(ij,2) = ij -1
        end if
        if (ij <=L) then
            NEIGHBORS(ij,3) = ij + (L-1)*L
        else
            NEIGHBORS(ij,3) =ij - L
        end if

        if (ij >= L*(L-1)) then
            NEIGHBORS(ij,4) = ij - (L-1)*L
        else
            NEIGHBORS(ij,4) = ij +L
        end if

      end do 
    end subroutine initialize_neighbors

    subroutine initialize_walls(LEFT_WALL, RIGHT_WALL, L)
        implicit none
        integer, intent(out) :: LEFT_WALL(:), RIGHT_WALL(:)
        integer, intent(in) :: L
        integer :: i
        LEFT_WALL(1) = 1
        RIGHT_WALL(1) = L

        do i = 2, L
            LEFT_WALL(i) = LEFT_WALL(i-1) + L
            RIGHT_WALL(i) = RIGHT_WALL(i-1) + L
        end do

    end subroutine initialize_walls

    subroutine initialize_lattice(LATTICE_SITES, mute_indices2, N, L)
        implicit none
        integer, intent(out) :: LATTICE_SITES(:), mute_indices2(:)
        integer, intent(in) :: N, L
        integer :: i, j_rng, temp
        real :: rand_val
        LATTICE_SITES = 0

       do i = 1, L**2
            mute_indices2(i) = i
        end do

        do i = L**2, 2, -1
            call random_number(rand_val)
            j_rng = int(rand_val * i) + 1
            temp = mute_indices2(i)
            mute_indices2(i) = mute_indices2(j_rng)
            mute_indices2(j_rng) = temp
        end do
    end subroutine initialize_lattice

    subroutine initialize_velocities(v0s, N)
        implicit none
        integer, intent(out) :: v0s(:)
        integer, intent(in) :: N
        integer :: i
        real :: rand_val

        do i = 1, N
            call random_number(rand_val)
            v0s(i) = int(rand_val * 4) + 1
        end do

    end subroutine initialize_velocities

    subroutine place_particles(LATTICE_SITES, OCCUPIED_SITES, mute_indices2, N)
        implicit none
        integer, intent(inout) :: LATTICE_SITES(:), OCCUPIED_SITES(:)
        integer, intent(in) :: mute_indices2(:), N
        integer :: i, delta

        do i = 1, N
            delta = mute_indices2(i)
            OCCUPIED_SITES(i) = delta
            LATTICE_SITES(delta) = 1
        end do

    end subroutine place_particles

    subroutine shuffle_indices(mute_indices, N)
        implicit none
        integer, intent(inout) :: mute_indices(:)
        integer, intent(in) :: N
        integer :: i, j_rng, temp
        real :: rand_val

        do i = N, 2, -1
            call random_number(rand_val)
            j_rng = int(rand_val * i) + 1
            temp = mute_indices(i)
            mute_indices(i) = mute_indices(j_rng)
            mute_indices(j_rng) = temp
        end do
    end subroutine shuffle_indices


    subroutine get_klist(NEIGHBORS,OCCUPIED_SITES,NLIST,N,LATTICE_SITES)
        implicit none
        integer, intent(in) :: N
        integer, intent(in) :: LATTICE_SITES(:)
        integer, intent(in) :: NEIGHBORS(:,:)

        integer, intent(in) :: OCCUPIED_SITES(:)
        intrinsic :: findloc
        integer, intent(inout) :: NLIST(:,:)
        integer :: ik 
        integer :: n_idxs(4)
        integer :: site_now
        integer :: ni
        integer :: idnow
        integer :: nij
        do ik = 1,N

            site_now = OCCUPIED_SITES(ik)
            n_idxs = NEIGHBORS(site_now,:)
            do ni = 1,4

                nij = n_idxs(ni)
                idnow  = findloc(OCCUPIED_SITES, value= nij,dim=1) 

                if (LATTICE_SITES(nij)==1) then
                    NLIST(ik,ni) = idnow
                else
                    NLIST(ik,ni) = 0
                end if
            end do
        end do

    end subroutine get_klist



    subroutine get_cluster(NLIST,N,CLIST)


        implicit none
        integer,intent(in) :: N
        integer, intent(inout) :: NLIST(:,:)
        integer, intent(inout) :: CLIST(:)
        integer :: ii, jj, cluster_id, site
        logical :: visited(N)
        integer :: queue(N), head, tail
        integer :: clusters(N)

        visited = .false.
        clusters = -1
        cluster_id = 0

        do ii = 1, N
            if (.not. visited(ii)) then
              cluster_id = cluster_id + 1
              head = 1
              tail = 1
             
              queue(head) = ii
              visited(ii) = .true.
              CLIST(ii) = cluster_id
              
              do while (head <= tail)
                site = queue(head)
                head = head + 1
                do jj = 1, 4
                  
                  if (NLIST(site, jj) /= 0 .and. .not. visited(NLIST(site, jj))) then
                    !print *, NEIGHBORS(site, jj) 
                    tail = tail + 1
                  
                    queue(tail) =NLIST(site, jj)
               
                    visited(NLIST(site, jj)) = .true.
                    CLIST(NLIST(site, jj)) = cluster_id
                  end if
                end do
              end do
            end if
          end do
        

    end subroutine get_cluster
    
    subroutine collect_unique_values(CLIST, N, UNIQUE_LIST, num_unique, dummy_p,dummy_s)
        implicit none
        integer, intent(in) :: N
        integer, intent(in) :: CLIST(N)
        integer, intent(out) :: UNIQUE_LIST(N)
        integer, intent(out) :: num_unique
        integer :: COUNT_LIST(N)
        integer :: i_clist, j_clist
        logical :: is_unique
        real(8), intent(inout) :: dummy_p
        real(8), intent(inout) :: dummy_s
        integer :: count_sum
        num_unique = 0
        COUNT_LIST = 0
        dummy_p = 0
        dummy_s = 0
        do i_clist = 1, N
            is_unique = .true.
           
            do j_clist = 1, num_unique
                if (CLIST(i_clist) == UNIQUE_LIST(j_clist)) then
                    is_unique = .false.
                    COUNT_LIST(j_clist) = COUNT_LIST(j_clist) + 1
                    exit
                end if
            end do
         
            if (is_unique) then
                num_unique = num_unique + 1
                UNIQUE_LIST(num_unique) = CLIST(i_clist)
                COUNT_LIST(num_unique) = 1
            end if
        end do
        do i_clist = 1, num_unique

            if (COUNT_LIST(i_clist) > dummy_p) then
                dummy_p = COUNT_LIST(i_clist)
            end if
            count_sum = count_sum + COUNT_LIST(i_clist)
        end do
        if (num_unique > 0) then
            dummy_s = real(count_sum) / real(num_unique)
        end if

    
    end subroutine collect_unique_values


end module simulation_utils




module v0_module
    implicit none
contains
subroutine update_pe(alpha,current_dir,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
  
    real,intent(in) :: alpha
    integer :: temp_int
   
    integer, intent(inout) :: current_dir
    real :: rand_val
    integer :: rand_int2
    integer,intent(in) :: L 
    integer :: old_delta
    integer, dimension(:),intent(in) :: RIGHT_WALL
    integer,dimension(:),intent(in)  :: LEFT_WALL
    integer, dimension(:), intent(inout) :: LATTICE_SITES
    integer, dimension(:), intent(inout) :: OCCUPIED_SITES
    integer, intent(in) :: idx
    integer :: new_delta_pt1
    integer :: new_delta_pt2
    integer :: new_delta
    integer :: Delta1, Delta2, Delta3, Delta4

    call random_number(rand_val)
  
    if (rand_val < alpha) then

        do 
        call random_number(rand_val)
        temp_int = int(rand_val * (4)) + 1
            if (temp_int /= current_dir) exit
        end do 

        current_dir = temp_int
    end if

    old_delta = OCCUPIED_SITES(idx)
   
    Delta1 = old_delta + (-L)*(1 - sign(1, L - old_delta))/2    + ( (L-1)*L)*( 1 + sign(1, L-old_delta))/2

    
    Delta2 = old_delta +1 - L *(1 - sign(1, mod(old_delta, L)-1))/2

    Delta3 = old_delta + L*(1 - sign(1, old_delta - ( L * (L-1) +1 )))/2   + (-L*(L-1))*(1+sign(1, old_delta - ( L * (L-1) +1 )))/2
    Delta4 = old_delta -1 + (L-1)*(1 - sign(1,mod(old_delta -1 ,L)-1))/2

    new_delta_pt1 = Delta1*(1 +merge(1,-1,current_dir==1))/2 + Delta2*(1 +merge(1,-1,current_dir==2))/2 
    new_delta_pt2= Delta3*(1 +merge(1,-1,current_dir==3))/2 + Delta4*(1 +merge(1,-1,current_dir==4))/2

    new_delta = new_delta_pt1 + new_delta_pt2
    
    LATTICE_SITES(old_delta) = LATTICE_SITES(old_delta) - (1 - LATTICE_SITES(new_delta))
    OCCUPIED_SITES(idx) = new_delta*(1 - LATTICE_SITES(new_delta)) + old_delta*LATTICE_SITES(old_delta)

    
    !old_delta = OCCUPIED_SITES(idx)
    LATTICE_SITES(new_delta) = LATTICE_SITES(new_delta) + (1 - LATTICE_SITES(old_delta))
  
!    old_delta = OCCUPIED_SITES(idx)


    return 
end subroutine update_pe
end module v0_module


module beta_module
    implicit none
contains
subroutine update_beta(beta,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
  
    real,intent(in) :: beta
    integer :: temp_int
   
    
    real :: rand_val
    integer :: rand_int2
    integer,intent(in) :: L 
    integer :: old_delta
    integer, dimension(:),intent(in) :: RIGHT_WALL
    integer,dimension(:),intent(in)  :: LEFT_WALL
    integer, dimension(:), intent(inout) :: LATTICE_SITES
    integer, dimension(:), intent(inout) :: OCCUPIED_SITES
    integer, intent(in) :: idx
    integer :: new_delta_pt1
    integer :: new_delta_pt2
    integer :: new_delta
    integer :: Delta1, Delta2, Delta3, Delta4

    call random_number(rand_val)
  
    
    temp_int = int(rand_val * (4)) + 1
  

    old_delta = OCCUPIED_SITES(idx)
   
    Delta1 = old_delta + (-L)*(1 - sign(1, L - old_delta))/2    + ( (L-1)*L)*( 1 + sign(1, L-old_delta))/2

    
    Delta2 = old_delta +1 - L *(1 - sign(1, mod(old_delta, L)-1))/2

    Delta3 = old_delta + L*(1 - sign(1, old_delta - ( L * (L-1) +1 )))/2   + (-L*(L-1))*(1+sign(1, old_delta - ( L * (L-1) +1 )))/2
    Delta4 = old_delta -1 + (L-1)*(1 - sign(1,mod(old_delta -1 ,L)-1))/2

    new_delta_pt1 = Delta1*(1 +merge(1,-1,temp_int==1))/2 + Delta2*(1 +merge(1,-1,temp_int==2))/2 
    new_delta_pt2= Delta3*(1 +merge(1,-1,temp_int==3))/2 + Delta4*(1 +merge(1,-1,temp_int==4))/2

    new_delta = new_delta_pt1 + new_delta_pt2
    
    LATTICE_SITES(old_delta) = LATTICE_SITES(old_delta) - (1 - LATTICE_SITES(new_delta))
    OCCUPIED_SITES(idx) = new_delta*(1 - LATTICE_SITES(new_delta)) + old_delta*LATTICE_SITES(old_delta)
  
    !old_delta = OCCUPIED_SITES(idx)
    LATTICE_SITES(new_delta) = LATTICE_SITES(new_delta) + (1 - LATTICE_SITES(old_delta))
  
!    old_delta = OCCUPIED_SITES(idx)


    return 
end subroutine update_beta
end module beta_module

program rnt_msd
    !guiblo = (gui)lherm + pa (blo) :)
  !  use traj_module
    use v0_module
    use beta_module
    use simulation_utils
    use omp_lib
    implicit none
    real, parameter :: alpha_values(1) = [0.001] ! valores baixos de alpha para comparar com o resultado contínuo
  !  real,parameter :: beta_values(50) = [0.0, 0.004,0.008 , &
   ! 0.012,0.015,0.018, 0.01863158,0.01926316, 0.01989474 ,  0.02052632,0.02115789 ,&
   ! 0.02178947,0.02242105,0.02305263, 0.02368421,0.02431579,0.02494737,0.02557895 ,&
   ! 0.02621053,0.02684211,0.02747368,0.02810526,0.02873684,0.02936842, 0.03,0.03291667,&
   !  0.03583333 ,0.03875 ,0.04166667,0.04458333,0.0475,0.05041667,0.05333333, 0.05625 ,&
   !  0.05916667,0.06208333,0.065,0.06791667,0.07083333,0.07375,0.07666667,0.07958333,0.0825, &    
   !  0.08541667, 0.08833333,0.09125,0.09416667,0.09708333,0.1,0.11]

    

    integer, parameter :: steps = 10000
    integer, parameter :: L = 96
  !  real :: phi_values() 
  !  phi_values(1)=0.04
 
    
    real :: phi
   ! real, parameter:: phi =0.1
    
    integer:: N
    
    !integer, parameter :: N = 1
    integer :: i 
    
    integer:: j
    integer :: k
    integer :: m
    real :: alpha ! run and tumble
     real :: beta ! noise
    integer :: step 
    real :: rand_val ! pro rng
    integer :: LATTICE_SITES(L**2) 
    integer ,allocatable :: OCCUPIED_SITES(:)
    integer :: RIGHT_WALL(L)
    integer:: NEIGHBORS(L**2,4)
    integer :: LEFT_WALL(L) 
    integer :: rand_int 
    integer :: rand_int2
    character(len=50) :: file_name 
    integer :: ni ! pra mexer no msd no fim do código 
    integer ,allocatable :: mute_indices(:)
    integer :: mute_indices2(L**2)
    integer :: temp
    integer :: MCSTEPS = 3000000
    integer,parameter :: COLSTEPS = 5000000
    integer :: i_now 
    integer :: j_rng
    integer :: i_rng
    integer :: delta
    integer :: old_delta
    integer :: idx
    integer,allocatable  :: v0s(:)
    integer :: current_dir
    integer :: noise_first
    integer ,allocatable :: NLIST(:,:)
    integer,allocatable  :: CLIST(:)
    integer,allocatable  :: UNIQUE_LIST(:)
    integer,parameter :: n_samples = 10000
    integer :: s
    integer :: count_times
    integer :: il
    integer :: num_unique
    integer, parameter :: num_betas = 60
    real(8) :: p 
    real(8) :: p2
    real(8) :: p4 
    real(8) :: dummy_p 
    real(8) :: dummy_s 
    real(8) :: av_s
    real(8) :: fmax(num_betas)
    real(8) :: fmax2(num_betas)
    real(8) :: fmax4(num_betas)
    real(8) :: u4(num_betas)
    real(8) :: meancl(num_betas)
    real(8) :: betas_print(num_betas)
    integer :: check_point = COLSTEPS/n_samples

   call random_seed()
   LEFT_WALL(1) = 1
   RIGHT_WALL(1) = L
   call initialize_walls(LEFT_WALL, RIGHT_WALL, L)
   call initialize_neighbors(LEFT_WALL,RIGHT_WALL,L,NEIGHBORS)
    
    phi =0.2
    N = int(phi * L * L)
    do i = 1,1
        alpha = 0.01
       
        
        allocate(OCCUPIED_SITES(N), mute_indices(N), v0s(N))
        allocate(NLIST(N, 4), CLIST(N),UNIQUE_LIST(N))
        beta = 0.0
        do j = 1,num_betas
	    print *, alpha
	    	    
            beta = beta + 0.005
            p =0
            p2  = 0 
            p4 = 0 
            av_s =0 
            dummy_s =0 
            dummy_p =0
            betas_print(j) = beta
            print *, beta
            call random_seed()
            call initialize_lattice(LATTICE_SITES, mute_indices2, N, L)
            call initialize_velocities(v0s, N)
            
            call place_particles(LATTICE_SITES, OCCUPIED_SITES, mute_indices2, N)

            do step = 1, MCSTEPS
                call random_seed()

                    do i_rng = 1, N
                        mute_indices(i_rng) = i_rng 
                    end do
                    
                    call shuffle_indices(mute_indices, N)
                  
                    do k = 1, N
                        noise_first = 0
                        idx = mute_indices(k)
                        current_dir = v0s(idx)

                        call random_number(rand_val)

                        if (rand_val < 0.5) then
                            call update_pe(alpha,current_dir,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            v0s(idx) = current_dir
                            call random_number(rand_val)
                            if (rand_val<beta) then
                            call update_beta(beta,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            end if 

                        else
                            call random_number(rand_val)
                            if (rand_val<beta) then
                            call update_beta(beta,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            end if 
                            call update_pe(alpha,current_dir,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            v0s(idx) = current_dir
                        end if 


                    end do       
            end do
	        count_times =0
            do step = 1, COLSTEPS
                call random_seed()
		    count_times = count_times +1
                    do i_rng = 1, N
                        mute_indices(i_rng) = i_rng 
                    end do
                    
                    call shuffle_indices(mute_indices, N)
                    
                    do k = 1, N
                        noise_first = 0
                        idx = mute_indices(k)
                       ! delta = OCCUPIED_SITES(idx)
                        current_dir = v0s(idx)
                      
                        call random_number(rand_val)

                        if (rand_val < 0.5) then
                            call update_pe(alpha,current_dir,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            v0s(idx) = current_dir
                            call random_number(rand_val)
                            if (rand_val<beta) then
                            call update_beta(beta,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            end if 

                        else
                            call random_number(rand_val)
                            if (rand_val<beta) then
                            call update_beta(beta,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            end if 
                            call update_pe(alpha,current_dir,RIGHT_WALL,LEFT_WALL,L,LATTICE_SITES,OCCUPIED_SITES,idx)
                            v0s(idx) = current_dir
                        end if 
                    end do    
if (count_times==check_point) then   
		    count_times=0
                    call get_klist(NEIGHBORS,OCCUPIED_SITES,NLIST,N,LATTICE_SITES)
                    call get_cluster(NLIST,N,CLIST)
                    call collect_unique_values(CLIST, N, UNIQUE_LIST, num_unique,dummy_p,dummy_s)

                   ! dummy_ p = dummy_p/N
                    av_s = av_s + dummy_s
                    p = p + dummy_p/N
                  !  print *, av_s
                 !   print *, p
                    p2 = p2 + (dummy_p/N)*(dummy_p/N)
                    p4 = p4 + (dummy_p/N)*(dummy_p/N)*(dummy_p/N)*(dummy_p/N)
 
			end if
            end do
       ! print *,p
     !  print *, av_s/n_samples
        fmax(j) = p/n_samples
        fmax2(j) = p2/n_samples
        fmax4(j) = p4/n_samples
        u4(j) =1  - (p4/n_samples)/(3*(p2/n_samples)*(p2/n_samples))
        meancl(j) = av_s/n_samples

       print *, fmax(j)
         

        end do
        deallocate(OCCUPIED_SITES, mute_indices, v0s, NLIST, CLIST)

    end do

    write(file_name, '(a,i0,a,i0,a,i0,a)') 'phi_', i, '_beta_',j, '_sample_',step,'.txt'

    open(unit=1, file=file_name)
       do ni = 1,num_betas

        write(1,*) fmax(ni),fmax2(ni),fmax4(ni), u4(ni), meancl(ni), betas_print(j)

    end do
    close(1)

   
end program rnt_msd

