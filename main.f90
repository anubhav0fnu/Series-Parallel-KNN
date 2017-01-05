Program To_Find_the_Knn

   use omp_lib
   implicit none 
   integer,parameter    :: n =500, m  = 11
   integer              :: class_1, class_2, TN, FP , FN , TP , no_class,nn
   real                 :: train(n,m), test(n,m), p1, p2, sensitivity, specificity
   integer              :: class_assgn(n),trainC(n), testC(n)
   integer              :: i, j, ktrial, prior,k,temp,p
   integer, parameter   :: ntrial=512
   double precision     :: flush(3*1024*1024) !flush caches
   double precision     :: t, t0, t1, time(0:ntrial), average
   
Setting the order of "K"
nn=int(n**0.6) !!optimazation used

!reading training set file and preparing train matrix
   open(unit=12,status='old',file='100032.dat') 
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
!reading test set file and preparing test matrix
   open(unit=17,status='old',file='c100032.dat')
   do i = 1,n
      read(17,2) testC(i),(test(i,j),j=1,m)
 2    format(i3,11f8.5)
   enddo

!calculation for K-nearest neighbours.
do k = 1,15

!calculating for both priors 
       p1 = float(class_1)/float(n)
       p2 = float(class_2)/float(n)
!calculation the time statistics   
  do ktrial=0,ntrial
  call random_number(flush)!flush caches
   t0 =omp_get_wtime()
   !***************************************
   call nearest_neighbours(train,test,n,m,k,class_assgn,trainC,p1,p2)
   !***************************************
   t1 =omp_get_wtime()
   time(ktrial)=(t1-t0)
  enddo
   average = sum(time(1:ntrial))/dble(ntrial)
 
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

   print*,'Results of     ',k,'- nearest neighbours'
   print*,'Priors: proportional,     p1= ',p1,'p2= ',p2
   print*,'----------------------------------------------------'
   print*,''  
   print*,'********************CLASSIFICATION******************'
   print*,''
   print*,'CONFUSION MATRIX->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
   print*,'----------------------------------------------------'
   print*,'                      |      PREDICTED CLASS        '
   print*,'                      | class1(-) | class2(+)| TOTAL'
   print*,'----------------------|---------- | ---------| -----'
   write(*,3) TN, FP, TN+FP
 3 FORMAT(' TRUE        CLASS1(-) |',i7,8x,i7,7x,i7)
   write(*,4) FN, TP, FN+TP
 4 FORMAT(' CLASS       CLASS2(+) |',i7,8x,i7,7x,i7)
   print*,'----------------------|---------- | ---------| -----'
   write(*,5) TN+FN, FP+TP, n
 5 format('                TOTAL  |',i7,8x,i7,7x,i7) 
   print*,'----------------------------------------------------'
   print*,''
   print*,'TOTAL NUMBER OF MIS-CLASSIFIED TEST OBSERVATION ',FP+FN
   print*,'TOTAL NUMBER OF CORRECTLY CLASSIFIED TEST OBSERVATION',TN+TP
   print*,''
   print*,'PERFORMANCE MEASURES OF CLASSIFICATION'
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
   print*,'NUMBER OF UNCLASSIFIED OBERVATIONS',no_class
   endif
   print*,'--------------EFFICIENCY OF CLASSIFICATION-----------------'
   print*,'PRECISION OF CLASSIFICATION',(float(TP)/float(TP+FP))
   print*,'ACCURACY OF CLASSIFICATION',(float(TP+TN)/float(n))
   print*,''
   print*,'AVERAGE TIME FOR SERIAL EXECUTION',average
enddo
end Program To_Find_the_KNN
