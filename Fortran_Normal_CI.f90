
!Felix Bittmann, 2019
!###############################################################################
real function fmean(arr, n) result(m)
	!Return the arithmetic mean of an array
	implicit none
	integer::n
	real, dimension(n)::arr
	m = sum(arr) / n
end function fmean
!###############################################################################
real function fstd(arr, n) result(m)
	!Returns the standard deviation of an array
	implicit none
	integer::n, i
	real::mean, fmean, summe
	real, dimension(n)::arr
	mean = fmean(arr, n)
	summe = 0
	do i = 1, n
		summe = summe + (arr(i) - mean)**2
	end do
	m = sqrt(summe / (n - 1))
end function fstd
!###############################################################################
real function fbootse(arr, n, reps) result(m)
	!Calculates the standard error of the mean for a given array
	implicit none
	integer::n, i, j, reps
	real, dimension(n)::sample, arr
	real, dimension(reps)::results
	real::fstd, r1, fmean
	iloop: do i = 1, reps
		jloop: do j = 1, n
			call random_number(r1)
			sample(j) = arr(int(r1 * n + 1))
		end do jloop
		results(i) = fmean(sample, n)			!Change desired function here
	end do iloop
	m = fstd(results, reps)
end function fbootse
!###############################################################################
SUBROUTINE init_random_seed2()
	!Call this function to reset the random number generator. Otherwise all results will be the same
	!https://stackoverflow.com/questions/31174367/slow-random-seed-generator-why
	!All Credit goes to the original author, Joel DeWitt
  USE ISO_Fortran_env, ONLY: INT64
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: seed(:)
  INTEGER :: i, n, un, istat, dt(8), pid
  INTEGER(INT64) :: t
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  OPEN(newunit=un, file='/dev/urandom', access='stream', status='old', action='read', form='unformatted', iostat=istat)
  IF (istat == 0) THEN
    READ(un) seed
    CLOSE(un)
  ELSE
    CALL SYSTEM_CLOCK(t)
    IF (t == 0) THEN
      CALL DATE_AND_TIME(values = dt)
      t = (dt(1) - 1970) * 365_INT64 * 24 * 60 * 60 * 1000 + dt(2) * 31_INT64 * 24 * 60 * 60 * 1000 + dt(3) * 24_INT64 * 60 * 60   &
      * 1000 + dt(5) * 60 * 60 * 1000 + dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)
    END IF
    pid = GETPID()
    t = IEOR(t, INT(pid, KIND(t)))
    DO i = 1, n
      seed(i) = lcg(t)
    END DO
  END IF
  CALL RANDOM_SEED(put = seed)
  DEALLOCATE(seed)
CONTAINS
  FUNCTION lcg(s)
    INTEGER :: lcg
    INTEGER(INT64) :: s
    IF (s == 0) THEN
      s = 104729
    ELSE
      s = MOD(s, 4294967296_INT64)
    END IF
    s = MOD(s * 279470273_INT64, 4294967291_INT64)
    lcg = INT(MOD(s, INT(HUGE(0), INT64)), KIND(0))
  END FUNCTION lcg
END SUBROUTINE init_random_seed2


PROGRAM  Main
	!Calculates a 95% normal based Confidence Interval for a given array
	IMPLICIT NONE
	real, DIMENSION(6) :: data1= (/1,2,3,4,5,6/)
	real::bootse, fbootse, fmean, mean
	call init_random_seed2()				!Call function to reset RNG
	bootse = fbootse(data1, 6, 50000)
	mean = fmean(data1, 6)
	write(*,*) mean - 1.96 * bootse, mean + 1.96 * bootse
	

END PROGRAM Main




