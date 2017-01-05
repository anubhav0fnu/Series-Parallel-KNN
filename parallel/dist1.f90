!Second Approach-(Using MPI + OpenMP)
subroutine distance_NN(train,test,n,m,test_row,k,send_class_NN,trainC,rank,p,chunk)
   use mpi
   implicit none
   integer                :: i, j, k, n, m, train_row, test_row, col, chunk
   real,intent(IN)        :: train(n,m), test(n,m)
   integer,intent(IN)     :: trainC(n)
   integer,intent(out)    :: send_class_NN(k)
   real                   :: dist(k),d,dist_1(k) 
   integer                :: class_NN(k),class_NN_1(k),p,rank,ierror,status(mpi_status_size)
   integer,parameter      :: rl = mpi_real,comm =mpi_comm_world
   integer                :: SEND_REQUEST(0:p-1,2),RECV_REQUEST(0:p-1,2)
   real                   :: recvbuf(chunk,m)
   
   call mpi_scatter(train,chunk*m,mpi_real,recvbuf,chunk*m,mpi_real,0, comm ,ierror)

!!!!!!!!!!!!!!!!!!!!!!calculations at rank 0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(rank==0)then
   !initialise the arrays
   call mpi_Irecv(send_class_NN(1),k,mpi_integer,p-1,12,comm,RECV_REQUEST,ierror)
   do i  = 1, k
      dist(i)     = 1.E+15
      class_NN(i) = 0
   enddo
!OpenMP implimentation:
!$omp parallel default(none) shared(dist,class_NN) private(train_row,i,j,col) firstprivate(test_row)     
!$omp do schedule(static)
   do train_row = 1,chunk
      d = 0.0
      do col = 1, m
	     !$omp  atomic update
         d = d + (recvbuf(train_row,col) - test(test_row,col))**2
      enddo
          do i = 1, k
         if (d <  dist(i)) then
            do j = k,i+1,-1
               if (j <= k) then
           !$omp  atomic update
                   dist(j)      = dist(j-1)
                   class_NN(j)  = class_NN(j-1)
               endif
            enddo
           !$omp  atomic update
         dist(i)      = d
         class_NN(i)  = trainC(train_row)
         goto 1
         endif
      enddo
 1       continue
   enddo
!$omp end do  
!$omp end parallel  
   call mpi_Isend(class_NN(1),k,mpi_integer,1,12,comm,ierror)
   CALL MPI_WAIT(SEND_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)
   call mpi_Isend(dist(1),k,rl,1,12,comm,SEND_REQUEST,ierror)
   CALL MPI_WAIT(SEND_REQUEST(rank,2),mpi_STATUS_ignore,IERROR)
   CALL MPI_WAIT(RECV_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)

!!!!!!!!!!!!!!!!!!!!!   calculations at last rank  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elseif(rank==p-1)then

   call mpi_Irecv(class_NN(1),k,mpi_integer,rank-1,12,comm,RECV_REQUEST,ierror)
   CALL MPI_WAIT(RECV_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)
   call mpi_Irecv(dist(1),k,rl,rank-1,12,comm,RECV_REQUEST,ierror)
   CALL MPI_WAIT(RECV_REQUEST(rank,2),mpi_STATUS_ignore,IERROR)
!OpenMP implimentation:
!$omp parallel default(none) shared(dist,class_NN) private(train_row,i,j,col) firstprivate(test_row)     
!$omp do schedule(static)
   do train_row = 1,chunk
      d = 0.0
      do col = 1, m
	  !$omp  atomic update
         d = d + (recvbuf(train_row,col) - test(test_row,col))**2
      enddo
          do i = 1, k
         if (d <  dist(i)) then
            do j = k,i+1,-1
               if (j <= k) then
			   !$omp  atomic update
                   dist(j)      = dist(j-1)
                   class_NN(j)  = class_NN(j-1)
               endif
            enddo
			!$omp  atomic update
         dist(i)      = d
         class_NN(i)  = trainC(train_row)
         goto 3
         endif
      enddo
 3       continue
   enddo
!$omp end do  
!$omp end parallel
   call mpi_Isend(class_NN(1),k,mpi_integer,0,12,comm,SEND_REQUEST,ierror)
   CALL MPI_WAIT(SEND_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)

!!!!!!!!!!!!!!!!!!!!!!!   calculations at all other ranks !!!!!!!!!!!!!!!!!!!!!!!!!!!
   else 
   call mpi_Irecv(class_NN_1(1),k,mpi_integer,rank-1,1,comm,RECV_REQUEST,ierror)
   CALL MPI_WAIT(RECV_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)
   call mpi_Irecv(dist_1(1),k,rl,rank-1,12,comm,RECV_REQUEST,ierror)
   CALL MPI_WAIT(RECV_REQUEST(rank,2),mpi_STATUS_ignore,IERROR)
!OpenMP implimentation:
!$omp parallel default(none) shared(dist,class_NN) private(train_row,i,j,col) firstprivate(test_row)     
!$omp do schedule(static) 
  do train_row = 1,chunk
      d = 0.0
      do col = 1, m
	  !$omp  atomic update
         d = d + (recvbuf(train_row,col) - test(test_row,col))**2
      enddo
          do i = 1, k
         if (d <  dist(i)) then
            do j = k,i+1,-1
               if (j <= k) then
			   !$omp  atomic update
                   dist_1(j)      = dist_1(j-1)
                   class_NN_1(j)  = class_NN_1(j-1)
               endif
            enddo
			!$omp  atomic update
         dist_1(i)      = d
         class_NN_1(i)  = trainC(train_row)
         goto 2
         endif
      enddo
 2       continue
   enddo
!$omp end do  
!$omp end parallel
   call mpi_Isend(class_NN_1(1),k,mpi_integer,rank+1,12,comm,SEND_REQUEST,ierror)
   CALL MPI_WAIT(SEND_REQUEST(rank,1),mpi_STATUS_ignore,IERROR)
   call mpi_Isend(dist_1(1),k,rl,rank+1,12,comm,SEND_REQUEST,ierror)
   CALL MPI_WAIT(SEND_REQUEST(rank,2),STATUS,IERROR)


   endif

end subroutine distance_NN
