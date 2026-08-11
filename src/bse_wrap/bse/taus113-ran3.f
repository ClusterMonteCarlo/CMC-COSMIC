C ====================================================================
C     This function is a drop-in replacement for the ran3 random
C     number generator from Numerical Recipes, Press et al.
C     It is based on E'cuyer's combined Tausworth generator with
C     period 2**113.
C
C ====================================================================
      real function ran3(IDUM)
        implicit none
        include "tausworth.h"
        integer(kind=int64) state(4), rnumber, taus113_gen_int
        integer i, IDUM, first, tid
c$      integer OMP_GET_THREAD_NUM
        COMMON/Taus113State/ state, first
        SAVE /Taus113State/
!$OMP THREADPRIVATE(/Taus113State/)
        EXTERNAL taus113_gen_int

        if(IDUM.lt.0) then
          call taus113_seeding(state, -IDUM)
C         The "warm-up" that has been left out in the seeding routine:
          print*, "Combined Tausworth 113 Generator warming up."
          do i=1,10
            call taus113_next_state(state)
          end do
          IDUM= -IDUM
        endif

        call taus113_next_state(state)
        rnumber= taus113_gen_int(state)

        if (first.eq.1) then
          tid= 1
c$        tid= OMP_GET_THREAD_NUM()
          first=0
          print("('taus-ran3: TID',I2,' The first random number is ',
     &    G16.8)"), tid, real(dble(rnumber)/dble(umax))
        endif

        ran3= real(dble(rnumber)/dble(umax))
      end function

      subroutine taus113_ran3_tester()
        integer i, failcount
        logical passed
        real rnumber, ran3
        external ran3

        passed=.true.
        failcount=0
        rnumber=ran3(-123456)
        do i=1,1000000000
          if (rnumber>=1.e0) then
            passed=.false.
            failcount= failcount + 1
            print*, "FOUL!!", rnumber-1.0e0, rnumber==1.e0, failcount, i
          end if
          rnumber=ran3(123456)
        end do
      end subroutine

C ====================================================================
C     RandomNormal and RandomTruncatedNormal — extracted from COSMIC's
C     src/cosmic/src/ran3.f. COSMIC's kick.f and assign_remnant.f call
C     these. We bundle them with our Tausworth ran3 so the internal
C     ran3() calls below resolve to the parallel-safe RNG (which is
C     critical correctness — NR's single-stream ran3 would correlate
C     draws across MPI ranks).
C ====================================================================

      SUBROUTINE RandomNormal(mean, sigma, idum, result)
* Generate a normally distributed random number with given mean and sigma
* using the Box-Muller transform

      real*8 mean, sigma, result
      integer idum
      real*8 u1, u2, Z0

      u1 = ran3(idum)
      u2 = ran3(idum)
      Z0 = SQRT(-2.d0*LOG(u1))*COS(2.d0*3.141592653589793d0*u2)
      result = Z0 * sigma + mean

      RETURN
      END


      SUBROUTINE RandomTruncatedNormal(mu, sigma, idum, lower, upper, x)
* Generate a random number from a truncated normal distribution
* with mean mu, standard deviation sigma, truncated to [lower, upper]
      IMPLICIT NONE
      REAL*8 mu, sigma, lower, upper, x
      INTEGER idum, max_attempts, attempt

      attempt = 0
      max_attempts = 1000

      do
          attempt = attempt + 1
          call RandomNormal(mu, sigma, idum, x)
          if (x .GE. lower .AND. x .LE. upper)then
             exit
          elseif (attempt.ge.max_attempts) then
            ! use the midpoint if we exceed max attempts
            x = 0.5d0 * (lower + upper)
            exit
          endif
      end do

      RETURN
      END
