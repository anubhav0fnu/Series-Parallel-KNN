Program To_Find_the_Knn_parallel
   use omp_lib
   use mpi
   implicit none 
   integer,parameter    :: n = 10, m  = 11
   integer              :: class_1, class_2, TN, FP , FN , TP , no_class,nn
   real                 :: train(n,m), test(n,m), p1, p2, sensitivity, specificity
   integer              :: class_assgn(n),trainC(n), testC(n)
   integer              :: i, j, ktrial, prior,k
   integer, parameter   :: ntrial=128
   double precision     :: flush(3*1024*1024) !flush caches
   double precision     :: t, t0, t1, time(0:ntrial),max_time(0:ntrial),average, gflops 
   integer,parameter    :: comm =  mpi_comm_world,dp=mpi_double_precision
   integer              :: ierror, p, rank

   call mpi_init(ierror)
   call mpi_comm_rank(comm,rank,ierror) 
   call mpi_comm_size(comm,p,ierror)
!Setting the order of "K"
nn=int(n**0.6) !!optimazation used
   
if(rank==0)then
!READING TRAINING FILE AND PREPARING TRAIN MATRIX
   open(unit=12,status='old',file='10.dat') 
   do i = 1, n
      read(12,1) trainC(i),(train(i,j),j=1,m)
 1    format(i3,11f8.5)
   enddo
!counting number of classes(1's and 2's) in training set.
   class_1 = 0
   class_2 = 0
   do i = 1,n
      if (trainC(i)==1)then 
         class_1 = class_1 + 1
      elseif(trainC(i)==2)then
         class_2 = class_2 + 1
      endif
   enddo
   print*,'*****************START*******************'
   print*,''
   print*,'Training set: cardinalities of 2 groups are',class_1, class_2
   print*,''
   print*,''
   print*,''
endif

!READING TEST FILE AND PREPARING TEST MATRIX ON ALL PROCESSORS
   open(unit=17,status='old',file='c10.dat')
   do i = 1,n
      read(17,2) testC(i),(test(i,j),j=1,m)
 2    format(i3,11f8.5)
   enddo

if(rank==0)then
k=nn
endif

calculating for both priors 

if(rank==0)then
   print*,''
  print*,'Results of',k,'- nearest neighbours'
   print*,'Priors: equal,       p1= ',p1,'p2= ',p2
   print*,'----------------------------------------------------'  
endif
!case of proportional priors
       p1 = float(class_1)/float(n)
       p2 = float(class_2)/float(n)
if(rank==0)then
   print*,''
   print*,'Results of     ',k,'- nearest neighbours'
   print*,'Priors: proportional,     p1= ',p1,'p2= ',p2
   print*,'----------------------------------------------------'
endif
!calculation the time statistics   
   do ktrial=0,ntrial
   call random_number(flush)!flush caches
   t0 =mpi_wtime()
   !***************************************
   call nearest_neighbours(train,test,n,m,k,class_assgn,trainC,p1,p2,rank,p)
   !***************************************
   t1 =mpi_wtime()
   time(ktrial)=(t1-t0)
   enddo
   call mpi_reduce(time,max_time,ntrial+1,dp,mpi_max,0,comm,ierror)

if(rank == 0)then
   print*,'Average time of Parallel version of KNN'
   average = sum(time(1:ntrial))/dble(ntrial)
  ! gflops = 2.d0*dble(n*n)*1.d-9/average
   TN = 0
   FP = 0
   FN = 0
   TP = 0
   no_class = 0
   do i = 1,n
    if     (class_assgn(i)==1 .and. testC(i)==1)then 
      TN = TN + 1
    elseif (class_assgn(i)==1 .and. testC(i)==2)then 
      FP = FP + 1
    elseif (class_assgn(i)==2 .and. testC(i)==1)then 
      FN = FN + 1
    elseif (class_assgn(i)==2 .and. testC(i)==2)then 
      TP = TP + 1
    elseif (class_assgn(i)==0)then 
      no_class = no_class + 1
    endif
   enddo
   print*,''  
   print*,'********************CLASSIFICATION******************'
   print*,''
   print*,'CONFUSION MATRIX->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
   print*,'----------------------------------------------------'
   print*,'                      |      PREDICTED CLASS        '
   print*,'                      | class1(-) | class2(+)| TOTAL'
   print*,'----------------------|---------- | ---------| -----'
   write(*,3) TN, FP, TN+FP
 3 FORMAT(' TRUE        CLASS1(-) |',i4,8x,i4,7x,i5)
   write(*,4) FN, TP, FN+TP
 4 FORMAT(' CLASS       CLASS2(+) |',i4,8x,i4,7x,i5)
   print*,'----------------------|---------- | ---------| -----'
   write(*,5) TN+FN, FP+TP, n
 5 format('                TOTAL  |',i4,8x,i4,7x,i5) 
   print*,'----------------------------------------------------'
   print*,''
   print*,'Important performance measures for classification'
   specificity=TN+FP
   if (specificity == 0)then
   print*,'SPECIFICITY =  cannot calculate'
   else
   print*,'SPECIFICITY = ',float(FP)/specificity
   endif
   sensitivity=FN+TP
   if (sensitivity == 0)then
   print*,'SENSITIVITY = cannot calculate'
   else
   print*,'SENSITIVITY = ',float(TP)/sensitivity
   endif
   if(no_class>0) then
   print*,'number of unclassified observations:',no_class
   endif
endif
  enddo--------------------------------------------------------------------->
if(rank==0)then
   print*,'Average time of PARALLEL execution',average
!   print*,'GFLOPS of PARALLEL execution',gflops
endif

call mpi_finalize(ierror)

end Program To_Find_the_KNN_parallel
