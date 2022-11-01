


subroutine geopack_recalc (iyear, idoy, ihour, imin, isec)    
    integer, intent(in) :: iyear, idoy, ihour, imin, isec
    call recalc (iyear, idoy, ihour, imin, isec)
end subroutine geopack_recalc


subroutine geopack_igrf_geo (r, theta, phi, br, btheta, bphi, ndata)
    integer :: ndata
    real, dimension(ndata), intent(in) :: r, theta, phi
    real, dimension(ndata), intent(out) :: br, btheta, bphi
    real :: bxout, byout, bzout, xin, yin, zin
    do n=1, ndata        
        rin = r(n)
        thetain = theta(n)
        phin = phi(n)
        call igrf_geo (rin, thetain, phin, brout, bthetaout, bphiout)
        br(n) = brout
        btheta(n) = bthetaout
        bphi(n) = bphiout
    end do
    
end subroutine geopack_igrf_geo


subroutine geopack_igrf_gsm (xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    real :: bxout, byout, bzout, xin, yin, zin
    
    do n=1,ndata
        
        xin = xgsm(n)
        yin = ygsm(n)
        zin = zgsm(n)
        call igrf_gsm (xin, yin, zin, bxout, byout, bzout)
        bxgsm(n) = bxout
        bygsm(n) = byout
        bzgsm(n) = bzout
        
    end do
end subroutine geopack_igrf_gsm


subroutine geopack_dip_gsm (xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    
    do n=1,ndata
        call dip (xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n))
    end do
end subroutine geopack_dip_gsm


subroutine geopack_sph2car (r, theta, phi, x, y, z, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: r, theta, phi    
    real, dimension(ndata), intent(out) :: x, y, z
    do n=1,ndata
        call sphcar (r(n), theta(n), phi(n), x(n), y(n), z(n), 1)
    end do
end subroutine geopack_sph2car


subroutine geopack_car2sph (r, theta, phi, x, y, z, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: r, theta, phi    
    real, dimension(ndata), intent(in) :: x, y, z
    do n=1,ndata
        call sphcar (r(n), theta(n), phi(n), x(n), y(n), z(n), -1)
    end do
end subroutine geopack_car2sph


subroutine geopack_bsph2car (theta, phi, br, btheta, bphi, bx, by, bz, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: theta, phi, br, btheta, bphi
    real, dimension(ndata), intent(out) :: bx, by, bz
    do n=1,ndata
        call bspcar (theta(n), phi(n), br(n), btheta(n), bphi(n), bx(n), by(n), bz(n))
    end do
end subroutine geopack_bsph2car


subroutine geopack_bcar2sph (x, y, z, bx, by, bz, br, btheta, bphi, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: x, y, z, bx, by, bz
    real, dimension(ndata), intent(out) :: br, btheta, bphi
    do n=1,ndata
        call bcarsp (x(n), y(n), z(n), bx(n), by(n), bz(n), br(n), btheta(n), bphi(n))
    end do
end subroutine geopack_bcar2sph


subroutine geopack_geo2mag (xgeo, ygeo, zgeo, xmag, ymag, zmag, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(out) :: xmag, ymag, zmag
    do n=1,ndata
        call geomag (xgeo(n), ygeo(n), zgeo(n), xmag(n), ymag(n), zmag(n), 1)
    end do
end subroutine geopack_geo2mag


subroutine geopack_mag2geo (xgeo, ygeo, zgeo, xmag, ymag, zmag, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(in) :: xmag, ymag, zmag
    do n=1,ndata
        call geomag (xgeo(n), ygeo(n), zgeo(n), xmag(n), ymag(n), zmag(n), -1)
    end do
end subroutine geopack_mag2geo


subroutine geopack_gei2geo (xgei, ygei, zgei, xgeo, ygeo, zgeo, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgei, ygei, zgei
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    do n=1,ndata
        call geigeo (xgei(n), ygei(n), zgei(n), xgeo(n), ygeo(n), zgeo(n), 1)
    end do
end subroutine geopack_gei2geo


subroutine geopack_geo2gei (xgei, ygei, zgei, xgeo, ygeo, zgeo, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgei, ygei, zgei
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    do n=1,ndata
        call geigeo (xgei(n), ygei(n), zgei(n), xgeo(n), ygeo(n), zgeo(n), -1)
    end do
end subroutine geopack_geo2gei


subroutine geopack_mag2sm (xmag, ymag, zmag, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xmag, ymag, zmag
    real, dimension(ndata), intent(out) :: xsm, ysm, zsm
    do n=1,ndata
        call magsm (xmag(n), ymag(n), zmag(n), xsm(n), ysm(n), zsm(n), 1)
    end do
end subroutine geopack_mag2sm


subroutine geopack_sm2mag (xmag, ymag, zmag, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xmag, ymag, zmag
    real, dimension(ndata), intent(in) :: xsm, ysm, zsm
    do n=1,ndata
        call magsm (xmag(n), ymag(n), zmag(n), xsm(n), ysm(n), zsm(n), -1)
    end do
end subroutine geopack_sm2mag


subroutine geopack_gsm2gse (xgsm, ygsm, zgsm, xgse, ygse, zgse, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: xgse, ygse, zgse
    do n=1,ndata
        call gsmgse (xgsm(n), ygsm(n), zgsm(n), xgse(n), ygse(n), zgse(n), 1)
    end do
end subroutine geopack_gsm2gse


subroutine geopack_gse2gsm (xgsm, ygsm, zgsm, xgse, ygse, zgse, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(in) :: xgse, ygse, zgse
    do n=1,ndata
        call gsmgse (xgsm(n), ygsm(n), zgsm(n), xgse(n), ygse(n), zgse(n), -1)
    end do
end subroutine geopack_gse2gsm


subroutine geopack_sm2gsm (xsm, ysm, zsm, xgsm, ygsm, zgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xsm, ysm, zsm
    real, dimension(ndata), intent(out) :: xgsm, ygsm, zgsm
    do n=1,ndata
        call smgsm (xsm(n), ysm(n), zsm(n), xgsm(n), ygsm(n), zgsm(n), 1)
    end do
end subroutine geopack_sm2gsm


subroutine geopack_gsm2sm (xsm, ysm, zsm, xgsm, ygsm, zgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xsm, ysm, zsm
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    do n=1,ndata
        call smgsm (xsm(n), ysm(n), zsm(n), xgsm(n), ygsm(n), zgsm(n), -1)
    end do
end subroutine geopack_gsm2sm


subroutine geopack_geo2gsm (xgeo, ygeo, zgeo, xgsm, ygsm, zgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(out) :: xgsm, ygsm, zgsm
    do n=1,ndata
        call geogsm (xgeo(n), ygeo(n), zgeo(n), xgsm(n), ygsm(n), zgsm(n), 1)
    end do
end subroutine geopack_geo2gsm


subroutine geopack_gsm2geo (xgeo, ygeo, zgeo, xgsm, ygsm, zgsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    do n=1,ndata
        call geogsm (xgeo(n), ygeo(n), zgeo(n), xgsm(n), ygsm(n), zgsm(n), -1)
    end do
end subroutine geopack_gsm2geo


subroutine geopack_trace (xi, yi, zi, dir, rlim, r0, iopt, parmod, exname, inname, xf, yf, zf, xfac, yfac, zfac, nfac)

    real, intent(in) :: xi, yi, zi
    real, intent(in) :: dir, rlim, r0, iopt
    real, dimension(10), intent(in) :: parmod
    character (len=*) :: exname,inname
    real, intent(out) :: xf, yf, zf
    real, dimension(1000),intent(out) :: xfac, yfac, zfac
    integer,intent(out) :: nfac
    
    call trace (xi, yi, zi, dir, rlim, r0, iopt, parmod, exname, inname, xf, yf, zf, xfac, yfac, zfac, nfac)

end subroutine geopack_trace


subroutine geopack_shuetal_mgnp (xn_pd, vel, bzimf, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, dist, id)

    real, intent(in) :: xn_pd, vel, bzimf, xgsm, ygsm, zgsm
    real, intent(out) :: xmgnp, ymgnp, zmgnp, dist
    integer, intent(out) :: id

    call shuetal_mgnp (xn_pd, vel, bzimf, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, dist, id)

end subroutine geopack_shuetal_mgnp


subroutine geopack_t96_mgnp (xn_pd, vel, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, dist, id)
    
    real, intent(in) :: xn_pd, vel, xgsm, ygsm, zgsm
    real, intent(out) :: xmgnp, ymgnp, zmgnp, dist
    integer, intent(out) :: id

    call t96_mgnp (xn_pd, vel, xgsm, ygsm, zgsm, xmgnp, ymgnp, zmgnp, dist, id)

end subroutine geopack_t96_mgnp



!    Now Geopack_2008 routines


subroutine geopack08_recalc (iyear, idoy, ihour, imin, isec, vgsex, vgsey, vgsez)
    integer, intent(in) :: iyear, idoy, ihour, imin, isec
    real, intent(in) :: vgsex, vgsey, vgsez
    call recalc_08 (iyear, idoy, ihour, imin, isec, vgsex, vgsey, vgsez)
end subroutine geopack08_recalc


subroutine geopack08_sun (iyear, idoy, ihour, imin, isec, gst, slong, srasn, sdec)
    integer, intent(in) :: iyear, idoy, ihour, imin, isec
    real, intent(out) :: gst, slong, srasn, sdec
    call sun_08 (iyear, idoy, ihour, imin, isec, gst, slong, srasn, sdec)
end subroutine geopack08_sun


subroutine geopack08_igrf_geo (r, theta, phi, br, btheta, bphi,ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: r, theta, phi
    real, dimension(ndata), intent(out) :: br, btheta, bphi
    do n=1,ndata
        call igrf_geo_08 (r(n), theta(n), phi(n), br(n), btheta(n), bphi(n))
    end do
end subroutine geopack08_igrf_geo


subroutine geopack08_igrf_gsw (xgsw, ygsw, zgsw, bxgsw, bygsw, bzgsw, ndata)

    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsw, ygsw, zgsw
    real, dimension(ndata), intent(out) :: bxgsw, bygsw, bzgsw
    
    do n=1,ndata
        call igrf_gsw_08 (xgsw(n), ygsw(n), zgsw(n), bxgsw(n), bygsw(n), bzgsw(n))
    end do
end subroutine geopack08_igrf_gsw


subroutine geopack08_dip_gsw (xgsw, ygsw, zgsw, bxgsw, bygsw, bzgsw, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsw, ygsw, zgsw
    real, dimension(ndata), intent(out) :: bxgsw, bygsw, bzgsw
    
    do n=1,ndata
        call dip_08 (xgsw(n), ygsw(n), zgsw(n), bxgsw(n), bygsw(n), bzgsw(n))
    end do
end subroutine geopack08_dip_gsw


subroutine geopack08_sph2car (r, theta, phi, x, y, z, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: r, theta, phi    
    real, dimension(ndata), intent(out) :: x, y, z
    do n=1,ndata
        call sphcar_08 (r(n), theta(n), phi(n), x(n), y(n), z(n), 1)
    end do
end subroutine geopack08_sph2car


subroutine geopack08_car2sph (r, theta, phi, x, y, z, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: r, theta, phi    
    real, dimension(ndata), intent(in) :: x, y, z
    do n=1,ndata
        call sphcar_08 (r(n), theta(n), phi(n), x(n), y(n), z(n), -1)
    end do
end subroutine geopack08_car2sph


subroutine geopack08_bsph2car (theta, phi, br, btheta, bphi, bx, by, bz, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: theta, phi, br, btheta, bphi
    real, dimension(ndata), intent(out) :: bx, by, bz
    do n=1,ndata
        call bspcar_08 (theta(n), phi(n), br(n), btheta(n), bphi(n), bx(n), by(n), bz(n))
    end do
end subroutine geopack08_bsph2car


subroutine geopack08_bcar2sph (x, y, z, bx, by, bz, br, btheta, bphi, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: x, y, z, bx, by, bz
    real, dimension(ndata), intent(out) :: br, btheta, bphi
    do n=1,ndata
        call bcarsp_08 (x(n), y(n), z(n), bx(n), by(n), bz(n), br(n), btheta(n), bphi(n))
    end do
end subroutine geopack08_bcar2sph


subroutine geopack08_gsw2gse (xgsw, ygsw, zgsw, xgse, ygse, zgse, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgsw, ygsw, zgsw
    real, dimension(ndata), intent(out) :: xgse, ygse, zgse
    do n=1,ndata
        call gswgse_08 (xgsw(n), ygsw(n), zgsw(n), xgse(n), ygse(n), zgse(n), 1)
    end do
end subroutine geopack08_gsw2gse


subroutine geopack08_gse2gsw (xgsw, ygsw, zgsw, xgse, ygse, zgse, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgsw, ygsw, zgsw
    real, dimension(ndata), intent(in) :: xgse, ygse, zgse
    do n=1,ndata
        call gswgse_08 (xgsw(n), ygsw(n), zgsw(n), xgse(n), ygse(n), zgse(n), -1)
    end do
end subroutine geopack08_gse2gsw


subroutine geopack08_geo2mag (xgeo, ygeo, zgeo, xmag, ymag, zmag, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(out) :: xmag, ymag, zmag
    do n=1,ndata
        call geomag_08 (xgeo(n), ygeo(n), zgeo(n), xmag(n), ymag(n), zmag(n), 1)
    end do
end subroutine geopack08_geo2mag


subroutine geopack08_mag2geo (xgeo, ygeo, zgeo, xmag, ymag, zmag, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(in) :: xmag, ymag, zmag
    do n=1,ndata
        call geomag_08 (xgeo(n), ygeo(n), zgeo(n), xmag(n), ymag(n), zmag(n), -1)
    end do
end subroutine geopack08_mag2geo


subroutine geopack08_gei2geo (xgei, ygei, zgei, xgeo, ygeo, zgeo, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgei, ygei, zgei
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    do n=1,ndata
        call geigeo_08 (xgei(n), ygei(n), zgei(n), xgeo(n), ygeo(n), zgeo(n), 1)
    end do
end subroutine geopack08_gei2geo


subroutine geopack08_geo2gei (xgei, ygei, zgei, xgeo, ygeo, zgeo, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgei, ygei, zgei
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    do n=1,ndata
        call geigeo_08 (xgei(n), ygei(n), zgei(n), xgeo(n), ygeo(n), zgeo(n), -1)
    end do
end subroutine geopack08_geo2gei


subroutine geopack08_mag2sm (xmag, ymag, zmag, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xmag, ymag, zmag
    real, dimension(ndata), intent(out) :: xsm, ysm, zsm
    do n=1,ndata
        call magsm_08 (xmag(n), ymag(n), zmag(n), xsm(n), ysm(n), zsm(n), 1)
    end do
end subroutine geopack08_mag2sm


subroutine geopack08_sm2mag (xmag, ymag, zmag, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xmag, ymag, zmag
    real, dimension(ndata), intent(in) :: xsm, ysm, zsm
    do n=1,ndata
        call magsm_08 (xmag(n), ymag(n), zmag(n), xsm(n), ysm(n), zsm(n), -1)
    end do
end subroutine geopack08_sm2mag


subroutine geopack08_sm2gsw (xgsw, ygsw, zgsw, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xsm, ysm, zsm
    real, dimension(ndata), intent(out) :: xgsw, ygsw, zgsw
    do n=1,ndata
        call smgsw_08 (xsm(n), ysm(n), zsm(n), xgsw(n), ygsw(n), zgsw(n),1)
    end do
end subroutine geopack08_sm2gsw


subroutine geopack08_gsw2sm (xgsw, ygsw, zgsw, xsm, ysm, zsm, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xsm, ysm, zsm
    real, dimension(ndata), intent(in) :: xgsw, ygsw, zgsw
    do n=1,ndata
        call smgsw_08 (xsm(n), ysm(n), zsm(n),xgsw(n), ygsw(n), zgsw(n),  -1)
    end do
end subroutine geopack08_gsw2sm


subroutine geopack08_geo2gsw (xgeo, ygeo, zgeo, xgsw, ygsw, zgsw, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(out) :: xgsw, ygsw, zgsw
    do n=1,ndata
        call geogsw_08 (xgeo(n), ygeo(n), zgeo(n), xgsw(n), ygsw(n), zgsw(n), 1)
    end do
end subroutine geopack08_geo2gsw


subroutine geopack08_gsw2geo (xgeo, ygeo, zgeo, xgsw, ygsw, zgsw, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: xgeo, ygeo, zgeo
    real, dimension(ndata), intent(in) :: xgsw, ygsw, zgsw
    do n=1,ndata
        call geogsw_08 (xgeo(n), ygeo(n), zgeo(n), xgsw(n), ygsw(n), zgsw(n), -1)
    end do
end subroutine geopack08_gsw2geo


subroutine geopack08_geod2geo (h, xmu, r, theta, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(in) :: h, xmu
    real, dimension(ndata), intent(out) :: r, theta
    do n=1,ndata
        call geodgeo_08 (h(n), xmu(n), r(n), theta(n), 1)
    end do
end subroutine geopack08_geod2geo


subroutine geopack08_geo2geod (h, xmu, r, theta, ndata)
    integer :: ndata, n
    real, dimension(ndata), intent(out) :: h, xmu
    real, dimension(ndata), intent(in) :: r, theta
    do n=1,ndata
        call geodgeo_08 (h(n), xmu(n), r(n), theta(n), -1)
    end do
end subroutine geopack08_geo2geod


subroutine geopack08_trace (xgsm, ygsm, zgsm, dir, dsmax, err, rlim, r0, iopt, parmod, &
                exname, inname, xf, yf, zf, xfac, yfac, zfac, nfac, nfacmax)

    real, intent(in) :: xgsm, ygsm, zgsm, dir, dsmax, err, rlim, r0, iopt    
    real, dimension(10), intent(in) :: parmod    
    character (len=*) :: exname, inname
    integer, intent(in) :: nfacmax
    real, intent(out) :: xf, yf, zf
    real, dimension(nfacmax),intent(out) :: xfac, yfac, zfac
    integer,intent(out) :: nfac
    
    call trace_08 (xgsm, ygsm, zgsm, dir, dsmax, err, rlim, r0, iopt, parmod, exname, inname, &
           xf, yf, zf, xfac, yfac, zfac, nfac, nfacmax)

end subroutine geopack08_trace


subroutine geopack08_shuetal_mgnp (xn_pd, vel, bzimf, xgsw, ygsw, zgsw, &
            xmgnp, ymgnp, zmgnp, dist, id)
    real, intent(in) :: xn_pd, vel, bzimf, xgsw, ygsw, zgsw
    real, intent(out) :: xmgnp, ymgnp, zmgnp, dist
    integer, intent(out) :: id
    call shuetal_mgnp_08 (xn_pd, vel, bzimf, xgsw, ygsw, zgsw, xmgnp, ymgnp, zmgnp, dist, id)
end subroutine geopack08_shuetal_mgnp


subroutine geopack08_t96_mgnp (xn_pd, vel, xgsw, ygsw, zgsw, xmgnp, ymgnp, zmgnp, dist, id)
    real, intent(in) :: xn_pd, vel, xgsw, ygsw, zgsw
    real, intent(out) :: xmgnp, ymgnp, zmgnp, dist
    integer, intent(out) :: id
    call t96_mgnp_08 (xn_pd, vel, xgsw, ygsw, zgsw, xmgnp, ymgnp, zmgnp, dist, id)
end subroutine geopack08_t96_mgnp


subroutine geopack_t89 (iopt, xgsm, ygsm, zgsm, bxgsm, bygsm ,bzgsm, ndata)
    integer :: ndata, n
    integer, intent(in) :: iopt
    real, dimension(10) :: parmod
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    COMMON /GEOPACK1/ A(15),PSI,AA(19)

    parmod = 0.0

    do n=1, ndata
        call t89d_sp (iopt, parmod, PSI, xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n))
    end do

end subroutine geopack_t89


subroutine geopack_t96 (parmod, xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, ndata)
    integer :: ndata, n, iopt
    real, dimension(10), intent(in) :: parmod
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    COMMON /GEOPACK1/ A(15),PSI,AA(19)

    iopt = 3

    do n=1, ndata
        call t96_01 (iopt, parmod, PSI, xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n))
    end do

end subroutine geopack_t96


subroutine geopack_t01 (parmod, xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, ndata)
    ! The t01_01c.f code has so many comments at the same lines after the code and variable declarations. 
    !   This causes the compiler to fail.  I have re-formatted the comments for the code compiles.
    integer :: ndata, n, iopt
    real, dimension(10), intent(in) :: parmod
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    real :: xtmp, ytmp, ztmp, bxtmp, bytmp, bztmp
    COMMON /GEOPACK1/ A(15),PSI,AA(19)

    iopt = 3

    do n=1, ndata
!        write(*,'(A, I2, A)'), '------------------', n ,'------------------------------'        
        xtmp = xgsm(n)
        ytmp = ygsm(n)
        ztmp = zgsm(n)
!        call t01_01 (iopt, parmod, PSI, xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n))
        call t01_01 (iopt, parmod, PSI, xtmp, ytmp, ztmp, bxtmp, bytmp, bztmp)
        bxgsm(n) = bxtmp
        bygsm(n) = bytmp
        bzgsm(n) = bztmp

!        write(*,'(I2, 6F10.3)'), n, bxtmp, bytmp, bztmp, xtmp, ytmp, ztmp
!        write(*,'(A)'), '------------------------------------------------'
    end do

end subroutine geopack_t01


subroutine geopack_t04 (parmod, xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, ndata)
    ! The TS04C model uses real*8 precisions. 
    !  This causes the compiler to fail. 
    integer :: ndata, n, iopt
    real, dimension(10), intent(in) :: parmod
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    COMMON /GEOPACK1/ A(15),PSI,AA(19)

    iopt = 3

    do n=1, ndata
        call T04_s (iopt, parmod, PSI, xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n))
    end do

end subroutine geopack_t04



subroutine geopack_bfield (iopt, parmod, xgsm, ygsm, zgsm, bxgsm, bygsm, bzgsm, exname, inname, ndata)
    integer :: ndata, n, iopt
    real, dimension(10), intent(in) :: parmod
    real, dimension(ndata), intent(in) :: xgsm, ygsm, zgsm
    real, dimension(ndata), intent(out) :: bxgsm, bygsm, bzgsm
    character (len=*) :: exname,inname
    COMMON /GEOPACK1/ A(15),PSI,AA(19)

    do n=1, ndata
        call bfield (xgsm(n), ygsm(n), zgsm(n), bxgsm(n), bygsm(n), bzgsm(n), iopt, parmod, exname, inname)        
    end do

end subroutine geopack_bfield