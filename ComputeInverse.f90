! Matrix Inverse Function - computes the inverse of the passed in matrix if it exits
! if the inverse does not exist, then the program will terminate

! Arguments
! A		REAL (n,n) (input)
!		User passed in nxn matrix 
!
! Ainv	REAL (n,n) (output)
!		Result of the inverse matrix operation
!
! n		INTEGER (input)
!		Dimension of the matrix

! Sources:
! DGETRF: http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
! DGETRI: http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga56d9c860ce4ce42ded7f914fdb0683ff.html

SUBROUTINE MATRIX_INVERSE(A, Ainv, n)
     INTEGER,  intent(in) :: n ! Dimension of the matrix
     REAL, dimension(n,n), intent(in) :: A ! Input matrix A
     REAL, dimension(n,n), intent(out) :: Ainv ! Inverse matrix of A

     ! Variables for the LAPACK inverse function
     REAL, dimension(size(A,1)) :: work ! Work Array
     INTEGER, dimension(size(A,1)) :: ipiv !	pivot indices
     INTEGER  :: info ! determines if inversion was successful

     ! LAPACK procedures
     external DGETRF
     external DGETRI

     Ainv = A ! Save the input matrix

     ! DGETRF computes an LU factorization of a general M-by-N matrix A
     ! using partial pivoting with row interchanges.
     call DGETRF(n,n,Ainv,n,ipiv,info)

     ! if info is not zero then the LU factorization was unsuccessful
     IF (info /= 0) THEN
          stop 'The input matrix is singular'
     END IF

     ! DGETRI computes the inverse of a matrix using the LU factorization
     ! computed by DGETRF.
     call DGETRI(n, Ainv,n,ipiv,work,n,info)

     IF (info /= 0) THEN
          stop 'The matrix inversion operation failed'
     END IF
     RETURN ! Terminate the function
END SUBROUTINE MATRIX_INVERSE


program Compute_Inverse
INTEGER :: n
REAL, dimension(4,4) :: A
REAL, dimension(4,4) :: Ainv

n = 4

!A = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 , 13, 14, 15 , 16 /), shape(A))
A = reshape((/ 5,6,6,8,2,2,2,8,6,6,2,8,2,3,6,7/), shape(A))
!A = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(A))

Print *,"Matrix A:"
Print *,A 

call MATRIX_INVERSE(A, Ainv, n)

Print *,"Matrix Ainv:"
Print *, Ainv

stop
end
