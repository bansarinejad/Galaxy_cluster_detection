      implicit real*8(a-h, o-z)
      double precision rinteg, rlfi
      dimension dm(1000), cum_lf(1000)
      
      ! Define cosmological parameters
      omega_m = 0.3d0
      c = 300000.0d0
      h = 50.0d0

      ! K-correction for different galaxy types (E/S0/Sab) in the i band
      rk = 2.5d0

      ! Schechter M* in Vega for E/S0/Sab in i band with H0=50, a=-alpha
      xms = -23.0d0
      a = 0.7d0
      dabsm = 0.02d0
      xmbr = xms - 4.0d0
      xmft = xms + 14.0d0
      nm = (xmft - xmbr) / dabsm

      ! i_vega band magnitude limit i_vega < 20.85mag (i_AB = 21.25)
      rmlim = 20.85d0
      rmbright = 12.0d0

      ! Open file to write results
      open(10, file='int_lf.txt')
      
      ! Define redshift parameters
      nz = 150
      zs = 0.0d0
      dz = 0.01
      
      ! Loop over redshift values
      do iz = 1, nz
         z = zs + dfloat(iz) * dz
         rcom = c / h * rinteg(z, omega_m)
         dl = rcom * (1.0d0 + z)
         dm(iz) = 25.0d0 + 5.0d0 * dlog10(dl) + rk * z
      enddo

      ! Write the magnitude limit to file
      write(10, 101) rmlim
  101 format(' redshift  Cumulative LF to mlim_Vega=', f6.2)
      
      ! Loop over redshift values again for cumulative luminosity function calculation
      do iz = 1, nz
         z = zs + dfloat(iz) * dz
         cum_lf(iz) = 0.0d0
         
         ! Loop over magnitude values
         do im = 1, nm
            xm = xmbr + dfloat(im) * dabsm
            apm = xm + dm(iz)
            ! Uncomment the following line for debugging if needed
            ! if(z.lt.0.1d0) write(*, *) z, apm, xm, dm(iz)
            if(apm.lt.rmlim.and.apm.ge.rmbright) then
               cum_lf(iz) = cum_lf(iz) + rlfi(xm, xms, a, xmbr, xmft)
            endif
         enddo
         
         write(10, 100) z, cum_lf(iz)
  100 format(f5.2, e16.6)
      enddo
      
      close(100)
      stop
      end

C Function to calculate RLFI
      double precision FUNCTION RLFI(XM, XMS, A, XMBR, XMFT)
      implicit real*8(a-h, o-z)
      if(XM.gt.XMBR.and.XM.lt.XMFT) then
         RLFI = 0.92d0 * dEXP(-dEXP(-0.92d0 * (XM - XMS)) + 0.92d0 * (A - 1.0d0) * (XM - XMS)) * 0.02d0
      else
         RLFI = 0.0d0
      endif
      return
      end

C Function to calculate RINTEG
      double precision FUNCTION RINTEG(Z, OMEGA)
      implicit real*8(a-h, o-z)
      DX = 0.00001d0
      NXMAX = (1.0d0 - 1.0d0 / (1.0d0 + Z)) / DX + 0.000001d0
      XS = 1.0d0 / (1.0d0 + Z) - DX / 2.0d0
      XINT = 0.0d0
      do IX = 1, NXMAX
         X = XS + FLOAT(IX) * DX
         F = ((1.0d0 - OMEGA + OMEGA / X + (1.0d0 - OMEGA) * (X * X - 1.0d0))**(-0.5d0)) / X
         XINT = XINT + F * DX
      enddo
      RINTEG = XINT
      return
      end
