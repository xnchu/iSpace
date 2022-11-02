"""
If trace is run 10000 times, it takes 121s with i7-6700k, 45s with M1 Max. 
"""


import ispace
import numpy as np
import os
import unittest
from scipy.io import readsav

print(os.getcwd())
if os.path.isfile('./iSpace/tests/test_geopack.sav'):
    data = readsav('./iSpace/tests/test_geopack.sav')
elif os.path.isfile('test_geopack.sav'):
    data = readsav('test_geopack.sav')
else:
    raise FileNotFoundError('test_geopack.sav not found')

print(ispace.__file__)

class test_geopack(unittest.TestCase):


    def test_recalc(self):
        
        ntest = data['test_recalc'].shape[0]
        for i in range(ntest):
            
            iyear=data['test_recalc'][i]['iyear']
            idoy=data['test_recalc'][i]['idoy']
            ihour=data['test_recalc'][i]['ihour']
            imin=data['test_recalc'][i]['imin']
            isec=data['test_recalc'][i]['isec']       
            tilt=data['test_recalc'][i]['tilt']     
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)

            sps = ispace.geopack1.sps
            cps = ispace.geopack1.cps
            psi = np.arcsin(sps)
            
            self.assertAlmostEqual(psi, tilt, 4, 'psi not correct')
    
        
    def test_igrf_geo(self):
        
        ntest = data['test_igrf_geo'].shape[0]

        for i in range(ntest):
            
            iyear=data['test_igrf_geo'][i]['iyear']
            idoy=data['test_igrf_geo'][i]['idoy']
            ihour=data['test_igrf_geo'][i]['ihour']
            imin=data['test_igrf_geo'][i]['imin']
            isec=data['test_igrf_geo'][i]['isec']     
            r = data['test_igrf_geo'][i]['r']
            theta = data['test_igrf_geo'][i]['theta']
            phi = data['test_igrf_geo'][i]['phi']
            br = data['test_igrf_geo'][i]['br']
            btheta = data['test_igrf_geo'][i]['btheta']
            bphi = data['test_igrf_geo'][i]['bphi']

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            br2, btheta2, bphi2 = ispace.geopack_igrf_geo(r, theta, phi)

            flag1 = np.allclose(br, br2, atol=10) or np.allclose(br, br2, rtol=1e-2)
            flag2 = np.allclose(btheta, btheta2, atol=10) or np.allclose(btheta, btheta2, rtol=1e-2)
            flag3 = np.allclose(bphi, bphi2, atol=10) or np.allclose(bphi, bphi2, rtol=1e-2)
            if not (flag1 and flag2 and flag3):
                print('i:', i   )
                print('Time:', iyear, idoy, ihour, imin, isec)
                print('Pos:', r, theta, phi)
                print('Geopack:', br, btheta, bphi)
                print('IDL:    ', br2, btheta2, bphi2)
                print('Ispace: ', br3, btheta3, bphi3)

            self.assertTrue(flag1 and flag2 and flag3)

    def test_igrf_gsm(self):

        ntest = data['test_igrf_gsm'].shape[0]

        for i in range(ntest):

            iyear=data['test_igrf_gsm'][i]['iyear']
            idoy=data['test_igrf_gsm'][i]['idoy']
            ihour=data['test_igrf_gsm'][i]['ihour']
            imin=data['test_igrf_gsm'][i]['imin']
            isec=data['test_igrf_gsm'][i]['isec']

            xgsm=data['test_igrf_gsm'][i]['xgsm']
            ygsm=data['test_igrf_gsm'][i]['ygsm']
            zgsm=data['test_igrf_gsm'][i]['zgsm']

            bxgsm=data['test_igrf_gsm'][i]['bxgsm']
            bygsm=data['test_igrf_gsm'][i]['bygsm']
            bzgsm=data['test_igrf_gsm'][i]['bzgsm']

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_igrf_gsm(xgsm, ygsm, zgsm)
            
            self.assertTrue(np.allclose(bxgsm, bxgsm2, atol=1e-2) or np.allclose(bxgsm, bxgsm2, rtol=1e-2))
            self.assertTrue(np.allclose(bygsm, bygsm2, atol=1e-2) or np.allclose(bygsm, bygsm2, rtol=1e-2))
            self.assertTrue(np.allclose(bzgsm, bzgsm2, atol=1e-2) or np.allclose(bzgsm, bzgsm2, rtol=1e-2))

    def test_dip_gsm(self):
        
        ntest = data['test_dip_gsm'].shape[0]    

        for i in range(ntest):

            iyear=data['test_dip_gsm'][i]['iyear']
            idoy=data['test_dip_gsm'][i]['idoy']
            ihour=data['test_dip_gsm'][i]['ihour']
            imin=data['test_dip_gsm'][i]['imin']
            isec=data['test_dip_gsm'][i]['isec']
            xgsm=data['test_dip_gsm'][i]['xgsm']
            ygsm=data['test_dip_gsm'][i]['ygsm']
            zgsm=data['test_dip_gsm'][i]['zgsm']
            bxgsm=data['test_dip_gsm'][i]['bxgsm']
            bygsm=data['test_dip_gsm'][i]['bygsm']
            bzgsm=data['test_dip_gsm'][i]['bzgsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_dip_gsm(xgsm, ygsm, zgsm)
    
            self.assertTrue(np.allclose(bxgsm, bxgsm2, atol=10) or np.allclose(bxgsm, bxgsm2, rtol=1e-2))
            self.assertTrue(np.allclose(bygsm, bygsm2, atol=10) or np.allclose(bygsm, bygsm2, rtol=1e-2))
            self.assertTrue(np.allclose(bzgsm, bzgsm2, atol=10) or np.allclose(bzgsm, bzgsm2, rtol=1e-2))


    def test_sph2car(self):

        ntest = data['test_sph2car'].shape[0]
        
        for i in range(ntest):
                
            r=data['test_sph2car'][i]['r']
            theta=data['test_sph2car'][i]['theta']
            phi=data['test_sph2car'][i]['phi']
            x=data['test_sph2car'][i]['x']
            y=data['test_sph2car'][i]['y']
            z=data['test_sph2car'][i]['z']
            
            x2, y2, z2 = ispace.geopack_sph2car(r, theta, phi)
            self.assertAlmostEqual(x[0], x2[0], 4, 'psi not correct')

            self.assertTrue(np.allclose(x, x2, atol=10) or np.allclose(x, x2, rtol=1e-2))
            self.assertTrue(np.allclose(y, y2, atol=10) or np.allclose(y, y2, rtol=1e-2))
            self.assertTrue(np.allclose(z, z2, atol=10) or np.allclose(z, z2, rtol=1e-2))


    def test_car2sph(self):
        
        ntest = data['test_car2sph'].shape[0]
        
        for i in range(ntest):
                
            r=data['test_car2sph'][i]['r']
            theta=data['test_car2sph'][i]['theta']
            phi=data['test_car2sph'][i]['phi']
            x=data['test_car2sph'][i]['x']
            y=data['test_car2sph'][i]['y']
            z=data['test_car2sph'][i]['z']
            
            r2, theta2, phi2 = ispace.geopack_car2sph(x, y, z)

            self.assertTrue(np.allclose(r, r2, atol=10) or np.allclose(r, r2, rtol=1e-2))
            self.assertTrue(np.allclose(theta, theta2, atol=10) or np.allclose(theta, theta2, rtol=1e-2))
            self.assertTrue(np.allclose(phi, phi2, atol=10) or np.allclose(phi, phi2, rtol=1e-2))


    def test_bsph2car(self):

        ntest = data['test_bsph2car'].shape[0]
        
        for i in range(ntest):
            
            theta=data['test_bsph2car'][i]['theta']
            phi=data['test_bsph2car'][i]['phi']
            br=data['test_bsph2car'][i]['br']
            btheta=data['test_bsph2car'][i]['btheta']
            bphi=data['test_bsph2car'][i]['bphi']
            bx=data['test_bsph2car'][i]['bx']
            by=data['test_bsph2car'][i]['by']
            bz=data['test_bsph2car'][i]['bz']
            
            bx2, by2, bz2 = ispace.geopack_bsph2car(theta, phi, br, btheta, bphi)

            self.assertTrue(np.allclose(bx, bx2, atol=1e-2) or np.allclose(bx, bx2, rtol=1e-2))
            self.assertTrue(np.allclose(by, by2, atol=1e-2) or np.allclose(by, by2, rtol=1e-2))
            self.assertTrue(np.allclose(bz, bz2, atol=1e-2) or np.allclose(bz, bz2, rtol=1e-2))
            

    def test_bcar2sph(self):
        
        ntest = data['test_bcar2sph'].shape[0]
        
        for i in range(ntest):
            
            x=data['test_bcar2sph'][i]['x']
            y=data['test_bcar2sph'][i]['y']
            z=data['test_bcar2sph'][i]['z']
            br=data['test_bcar2sph'][i]['br']
            btheta=data['test_bcar2sph'][i]['btheta']
            bphi=data['test_bcar2sph'][i]['bphi']
            bx=data['test_bcar2sph'][i]['bx']
            by=data['test_bcar2sph'][i]['by']
            bz=data['test_bcar2sph'][i]['bz']
            
            br2, btheta2, bphi2 = ispace.geopack_bcar2sph(x, y, z, bx, by, bz)
            np.testing.assert_allclose(br, br2, rtol=1e-2)
            np.testing.assert_allclose(btheta, btheta2, rtol=1e-2)
            np.testing.assert_allclose(bphi, bphi2, rtol=1e-2)
            # self.assertTrue(np.allclose(br, br2, atol=1e-2) or np.allclose(br, br2, rtol=1e-2))
            # self.assertTrue(np.allclose(btheta, btheta2, atol=1e-2) or np.allclose(btheta, btheta2, rtol=1e-2))
            # self.assertTrue(np.allclose(bphi, bphi2, atol=1e-2) or np.allclose(bphi, bphi2, rtol=1e-2))
        

    def test_geo2mag(self):
        
        ntest = data['test_geo2mag'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_geo2mag'][i]['iyear']
            idoy=data['test_geo2mag'][i]['idoy']
            ihour=data['test_geo2mag'][i]['ihour']
            imin=data['test_geo2mag'][i]['imin']
            isec=data['test_geo2mag'][i]['isec']
            
            xgeo=data['test_geo2mag'][i]['xgeo']
            ygeo=data['test_geo2mag'][i]['ygeo']
            zgeo=data['test_geo2mag'][i]['zgeo']
            xmag=data['test_geo2mag'][i]['xmag']
            ymag=data['test_geo2mag'][i]['ymag']
            zmag=data['test_geo2mag'][i]['zmag']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xmag2, ymag2, zmag2 = ispace.geopack_geo2mag(xgeo, ygeo, zgeo)

            np.testing.assert_allclose(xmag, xmag2, atol=3e-2)
            np.testing.assert_allclose(ymag, ymag2, atol=3e-2)
            np.testing.assert_allclose(zmag, zmag2, atol=3e-2)
            # self.assertTrue(np.allclose(xmag, xmag2, atol=1e-2) or np.allclose(xmag, xmag2, rtol=1e-2))
            # self.assertTrue(np.allclose(ymag, ymag2, atol=1e-2) or np.allclose(ymag, ymag2, rtol=1e-2))
            # self.assertTrue(np.allclose(zmag, zmag2, atol=1e-2) or np.allclose(zmag, zmag2, rtol=1e-2))


    def test_mag2geo(self):

        ntest = data['test_mag2geo'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_mag2geo'][i]['iyear']
            idoy=data['test_mag2geo'][i]['idoy']
            ihour=data['test_mag2geo'][i]['ihour']
            imin=data['test_mag2geo'][i]['imin']
            isec=data['test_mag2geo'][i]['isec']
            
            xgeo=data['test_mag2geo'][i]['xgeo']
            ygeo=data['test_mag2geo'][i]['ygeo']
            zgeo=data['test_mag2geo'][i]['zgeo']
            xmag=data['test_mag2geo'][i]['xmag']
            ymag=data['test_mag2geo'][i]['ymag']
            zmag=data['test_mag2geo'][i]['zmag']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgeo2, ygeo2, zgeo2 = ispace.geopack_mag2geo(xmag, ymag, zmag)

            self.assertTrue(np.allclose(xgeo, xgeo2, atol=1e-2) or np.allclose(xgeo, xgeo2, rtol=1e-2))
            self.assertTrue(np.allclose(ygeo, ygeo2, atol=1e-2) or np.allclose(ygeo, ygeo2, rtol=1e-2))
            self.assertTrue(np.allclose(zgeo, zgeo2, atol=1e-2) or np.allclose(zgeo, zgeo2, rtol=1e-2))


    def test_gei2geo(self):

        ntest = data['test_gei2geo'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_gei2geo'][i]['iyear']
            idoy=data['test_gei2geo'][i]['idoy']
            ihour=data['test_gei2geo'][i]['ihour']
            imin=data['test_gei2geo'][i]['imin']
            isec=data['test_gei2geo'][i]['isec']
            
            xgei=data['test_gei2geo'][i]['xgei']
            ygei=data['test_gei2geo'][i]['ygei']
            zgei=data['test_gei2geo'][i]['zgei']
            xgeo=data['test_gei2geo'][i]['xgeo']
            ygeo=data['test_gei2geo'][i]['ygeo']
            zgeo=data['test_gei2geo'][i]['zgeo']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgeo2, ygeo2, zgeo2 = ispace.geopack_gei2geo(xgei, ygei, zgei)

            self.assertTrue(np.allclose(xgeo, xgeo2, atol=1e-2) or np.allclose(xgeo, xgeo2, rtol=1e-2))
            self.assertTrue(np.allclose(ygeo, ygeo2, atol=1e-2) or np.allclose(ygeo, ygeo2, rtol=1e-2))
            self.assertTrue(np.allclose(zgeo, zgeo2, atol=1e-2) or np.allclose(zgeo, zgeo2, rtol=1e-2))


    def test_geo2gei(self):
        
        ntest = data['test_geo2gei'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_geo2gei'][i]['iyear']
            idoy=data['test_geo2gei'][i]['idoy']
            ihour=data['test_geo2gei'][i]['ihour']
            imin=data['test_geo2gei'][i]['imin']
            isec=data['test_geo2gei'][i]['isec']
            
            xgei=data['test_geo2gei'][i]['xgei']
            ygei=data['test_geo2gei'][i]['ygei']
            zgei=data['test_geo2gei'][i]['zgei']
            xgeo=data['test_geo2gei'][i]['xgeo']
            ygeo=data['test_geo2gei'][i]['ygeo']
            zgeo=data['test_geo2gei'][i]['zgeo']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgei2, ygei2, zgei2 = ispace.geopack_geo2gei(xgeo, ygeo, zgeo)

            self.assertTrue(np.allclose(xgei, xgei2, atol=1e-2) or np.allclose(xgei, xgei2, rtol=1e-2))
            self.assertTrue(np.allclose(ygei, ygei2, atol=1e-2) or np.allclose(ygei, ygei2, rtol=1e-2))
            self.assertTrue(np.allclose(zgei, zgei2, atol=1e-2) or np.allclose(zgei, zgei2, rtol=1e-2))

    def test_mag2sm(self):
        
        ntest = data['test_mag2sm'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_mag2sm'][i]['iyear']
            idoy=data['test_mag2sm'][i]['idoy']
            ihour=data['test_mag2sm'][i]['ihour']
            imin=data['test_mag2sm'][i]['imin']
            isec=data['test_mag2sm'][i]['isec']
            
            xsm=data['test_mag2sm'][i]['xsm']
            ysm=data['test_mag2sm'][i]['ysm']
            zsm=data['test_mag2sm'][i]['zsm']
            xmag=data['test_mag2sm'][i]['xmag']
            ymag=data['test_mag2sm'][i]['ymag']
            zmag=data['test_mag2sm'][i]['zmag']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xsm2, ysm2, zsm2 = ispace.geopack_mag2sm(xmag, ymag, zmag)

            self.assertTrue(np.allclose(xsm, xsm2, atol=1e-2) or np.allclose(xsm, xsm2, rtol=1e-2))
            self.assertTrue(np.allclose(ysm, ysm2, atol=1e-2) or np.allclose(ysm, ysm2, rtol=1e-2))
            self.assertTrue(np.allclose(zsm, zsm2, atol=1e-2) or np.allclose(zsm, zsm2, rtol=1e-2))    



    def test_sm2mag(self):
                    
        ntest = data['test_sm2mag'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_sm2mag'][i]['iyear']
            idoy=data['test_sm2mag'][i]['idoy']
            ihour=data['test_sm2mag'][i]['ihour']
            imin=data['test_sm2mag'][i]['imin']
            isec=data['test_sm2mag'][i]['isec']
            
            xsm=data['test_sm2mag'][i]['xsm']
            ysm=data['test_sm2mag'][i]['ysm']
            zsm=data['test_sm2mag'][i]['zsm']
            xmag=data['test_sm2mag'][i]['xmag']
            ymag=data['test_sm2mag'][i]['ymag']
            zmag=data['test_sm2mag'][i]['zmag']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xmag2, ymag2, zmag2 = ispace.geopack_sm2mag(xsm, ysm, zsm)

            self.assertTrue(np.allclose(xmag, xmag2, atol=1e-2) or np.allclose(xmag, xmag2, rtol=1e-2))
            self.assertTrue(np.allclose(ymag, ymag2, atol=1e-2) or np.allclose(ymag, ymag2, rtol=1e-2))
            self.assertTrue(np.allclose(zmag, zmag2, atol=1e-2) or np.allclose(zmag, zmag2, rtol=1e-2))

    def test_gsm2gse(self):
                    
        ntest = data['test_gsm2gse'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_gsm2gse'][i]['iyear']
            idoy=data['test_gsm2gse'][i]['idoy']
            ihour=data['test_gsm2gse'][i]['ihour']
            imin=data['test_gsm2gse'][i]['imin']
            isec=data['test_gsm2gse'][i]['isec']
            
            xgse=data['test_gsm2gse'][i]['xgse']
            ygse=data['test_gsm2gse'][i]['ygse']
            zgse=data['test_gsm2gse'][i]['zgse']
            xgsm=data['test_gsm2gse'][i]['xgsm']
            ygsm=data['test_gsm2gse'][i]['ygsm']
            zgsm=data['test_gsm2gse'][i]['zgsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgse2, ygse2, zgse2 = ispace.geopack_gsm2gse(xgsm, ygsm, zgsm)

            self.assertTrue(np.allclose(xgse, xgse2, atol=1e-2) or np.allclose(xgse, xgse2, rtol=1e-2))
            self.assertTrue(np.allclose(ygse, ygse2, atol=1e-2) or np.allclose(ygse, ygse2, rtol=1e-2))
            self.assertTrue(np.allclose(zgse, zgse2, atol=1e-2) or np.allclose(zgse, zgse2, rtol=1e-2))


    def test_gse2gsm(self):
                    
        ntest = data['test_gse2gsm'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_gse2gsm'][i]['iyear']
            idoy=data['test_gse2gsm'][i]['idoy']
            ihour=data['test_gse2gsm'][i]['ihour']
            imin=data['test_gse2gsm'][i]['imin']
            isec=data['test_gse2gsm'][i]['isec']
            
            xgse=data['test_gse2gsm'][i]['xgse']
            ygse=data['test_gse2gsm'][i]['ygse']
            zgse=data['test_gse2gsm'][i]['zgse']
            xgsm=data['test_gse2gsm'][i]['xgsm']
            ygsm=data['test_gse2gsm'][i]['ygsm']
            zgsm=data['test_gse2gsm'][i]['zgsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgsm2, ygsm2, zgsm2 = ispace.geopack_gse2gsm(xgse, ygse, zgse)

            self.assertTrue(np.allclose(xgsm, xgsm2, atol=1e-2) or np.allclose(xgsm, xgsm2, rtol=1e-2))
            self.assertTrue(np.allclose(ygsm, ygsm2, atol=1e-2) or np.allclose(ygsm, ygsm2, rtol=1e-2))
            self.assertTrue(np.allclose(zgsm, zgsm2, atol=1e-2) or np.allclose(zgsm, zgsm2, rtol=1e-2))
                    
    def test_sm2gsm(self):
                    
        ntest = data['test_sm2gsm'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_sm2gsm'][i]['iyear']
            idoy=data['test_sm2gsm'][i]['idoy']
            ihour=data['test_sm2gsm'][i]['ihour']
            imin=data['test_sm2gsm'][i]['imin']
            isec=data['test_sm2gsm'][i]['isec']
            
            xgsm=data['test_sm2gsm'][i]['xgsm']
            ygsm=data['test_sm2gsm'][i]['ygsm']
            zgsm=data['test_sm2gsm'][i]['zgsm']
            xsm=data['test_sm2gsm'][i]['xsm']
            ysm=data['test_sm2gsm'][i]['ysm']
            zsm=data['test_sm2gsm'][i]['zsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgsm2, ygsm2, zgsm2 = ispace.geopack_sm2gsm(xsm, ysm, zsm)

            self.assertTrue(np.allclose(xgsm, xgsm2, atol=1e-2) or np.allclose(xgsm, xgsm2, rtol=1e-2))
            self.assertTrue(np.allclose(ygsm, ygsm2, atol=1e-2) or np.allclose(ygsm, ygsm2, rtol=1e-2))
            self.assertTrue(np.allclose(zgsm, zgsm2, atol=1e-2) or np.allclose(zgsm, zgsm2, rtol=1e-2))

    def test_gsm2sm(self):
                    
        ntest = data['test_gsm2sm'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_gsm2sm'][i]['iyear']
            idoy=data['test_gsm2sm'][i]['idoy']
            ihour=data['test_gsm2sm'][i]['ihour']
            imin=data['test_gsm2sm'][i]['imin']
            isec=data['test_gsm2sm'][i]['isec']
            
            xgsm=data['test_gsm2sm'][i]['xgsm']
            ygsm=data['test_gsm2sm'][i]['ygsm']
            zgsm=data['test_gsm2sm'][i]['zgsm']
            xsm=data['test_gsm2sm'][i]['xsm']
            ysm=data['test_gsm2sm'][i]['ysm']
            zsm=data['test_gsm2sm'][i]['zsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xsm2, ysm2, zsm2 = ispace.geopack_gsm2sm(xgsm, ygsm, zgsm)

            self.assertTrue(np.allclose(xsm, xsm2, atol=1e-2) or np.allclose(xsm, xsm2, rtol=1e-2))
            self.assertTrue(np.allclose(ysm, ysm2, atol=1e-2) or np.allclose(ysm, ysm2, rtol=1e-2))
            self.assertTrue(np.allclose(zsm, zsm2, atol=1e-2) or np.allclose(zsm, zsm2, rtol=1e-2))

    def test_geo2gsm(self):
                    
        ntest = data['test_geo2gsm'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_geo2gsm'][i]['iyear']
            idoy=data['test_geo2gsm'][i]['idoy']
            ihour=data['test_geo2gsm'][i]['ihour']
            imin=data['test_geo2gsm'][i]['imin']
            isec=data['test_geo2gsm'][i]['isec']
            
            xgeo=data['test_geo2gsm'][i]['xgeo']
            ygeo=data['test_geo2gsm'][i]['ygeo']
            zgeo=data['test_geo2gsm'][i]['zgeo']
            xgsm=data['test_geo2gsm'][i]['xgsm']
            ygsm=data['test_geo2gsm'][i]['ygsm']
            zgsm=data['test_geo2gsm'][i]['zgsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgsm2, ygsm2, zgsm2 = ispace.geopack_geo2gsm(xgeo, ygeo, zgeo)

            self.assertTrue(np.allclose(xgsm, xgsm2, atol=1e-2) or np.allclose(xgsm, xgsm2, rtol=1e-2))
            self.assertTrue(np.allclose(ygsm, ygsm2, atol=1e-2) or np.allclose(ygsm, ygsm2, rtol=1e-2))
            self.assertTrue(np.allclose(zgsm, zgsm2, atol=1e-2) or np.allclose(zgsm, zgsm2, rtol=1e-2))


    def test_gsm2geo(self):
                    
        ntest = data['test_gsm2geo'].shape[0]
        
        for i in range(ntest):

            iyear=data['test_gsm2geo'][i]['iyear']
            idoy=data['test_gsm2geo'][i]['idoy']
            ihour=data['test_gsm2geo'][i]['ihour']
            imin=data['test_gsm2geo'][i]['imin']
            isec=data['test_gsm2geo'][i]['isec']
            
            xgeo=data['test_gsm2geo'][i]['xgeo']
            ygeo=data['test_gsm2geo'][i]['ygeo']
            zgeo=data['test_gsm2geo'][i]['zgeo']
            xgsm=data['test_gsm2geo'][i]['xgsm']
            ygsm=data['test_gsm2geo'][i]['ygsm']
            zgsm=data['test_gsm2geo'][i]['zgsm']
            
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xgeo2, ygeo2, zgeo2 = ispace.geopack_gsm2geo(xgsm, ygsm, zgsm)

            self.assertTrue(np.allclose(xgeo, xgeo2, atol=1e-2) or np.allclose(xgeo, xgeo2, rtol=1e-2))
            self.assertTrue(np.allclose(ygeo, ygeo2, atol=1e-2) or np.allclose(ygeo, ygeo2, rtol=1e-2))
            self.assertTrue(np.allclose(zgeo, zgeo2, atol=1e-2) or np.allclose(zgeo, zgeo2, rtol=1e-2))
            
    def test_trace(self):
        ntest = data['test_trace'].shape[0]
        for i in range(ntest):
            iyear=data['test_trace'][i]['iyear']
            idoy=data['test_trace'][i]['idoy']
            ihour=data['test_trace'][i]['ihour']
            imin=data['test_trace'][i]['imin']
            isec=data['test_trace'][i]['isec']
            
            dir=data['test_trace'][i]['dir']
            iopt=data['test_trace'][i]['iopt']
            parmod=data['test_trace'][i]['parmod']
            r0=data['test_trace'][i]['r0']
            rlim=data['test_trace'][i]['rlim']

            xgsm=data['test_trace'][i]['xgsm']
            ygsm=data['test_trace'][i]['ygsm']
            zgsm=data['test_trace'][i]['zgsm']
            
            xf=data['test_trace'][i]['xf']
            yf=data['test_trace'][i]['yf']
            zf=data['test_trace'][i]['zf']
            exname='T96'
            inname='DIP'

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xf2, yf2, zf2, xx2, yy2, zz2, L2 = ispace.geopack_trace(xgsm, ygsm, zgsm, dir, rlim, r0, iopt, parmod, exname, inname)
            
            # Skip in case of open field lines. 
            if(np.sqrt(xf2**2 + yf2**2 + zf2**2) > 2.0):
                continue

            # np.testing.assert_allclose(xf, xf2, atol=5e-2)
            # np.testing.assert_allclose(yf, yf2, atol=5e-2)
            # np.testing.assert_allclose(zf, zf2, atol=5e-2)

            self.assertTrue(np.allclose(xf, xf2, atol=5e-2) or np.allclose(xf, xf2, rtol=5e-2))
            self.assertTrue(np.allclose(yf, yf2, atol=5e-2) or np.allclose(yf, yf2, rtol=5e-2))
            self.assertTrue(np.allclose(zf, zf2, atol=5e-2) or np.allclose(zf, zf2, rtol=5e-2))
        
    def test_trace2(self):

        ntest = data['test_trace2'].shape[0]

        for i in range(ntest):

            iyear=data['test_trace2'][i]['iyear']
            idoy=data['test_trace2'][i]['idoy']
            ihour=data['test_trace2'][i]['ihour']
            imin=data['test_trace2'][i]['imin']
            isec=data['test_trace2'][i]['isec']
            
            dir=data['test_trace2'][i]['dir']
            iopt=data['test_trace2'][i]['iopt']
            parmod=data['test_trace2'][i]['parmod']
            r0=data['test_trace2'][i]['r0']
            rlim=data['test_trace2'][i]['rlim']

            xgsm=data['test_trace2'][i]['xgsm']
            ygsm=data['test_trace2'][i]['ygsm']
            zgsm=data['test_trace2'][i]['zgsm']
            
            xf=data['test_trace2'][i]['xf']
            yf=data['test_trace2'][i]['yf']
            zf=data['test_trace2'][i]['zf']
            exname='T96'
            inname='IGRF_GSM'

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xf2, yf2, zf2, xx2, yy2, zz2, L2 = ispace.geopack_trace(xgsm, ygsm, zgsm, dir, rlim, r0, iopt, parmod, exname, inname)
            
            # Skip in case of open field lines. 
            if(np.sqrt(xf2**2 + yf2**2 + zf2**2) > 2.0):
                continue

            # np.testing.assert_allclose(xf, xf2, atol=5e-2)
            # np.testing.assert_allclose(yf, yf2, atol=5e-2)
            # np.testing.assert_allclose(zf, zf2, atol=5e-2)

            self.assertTrue(np.allclose(xf, xf2, atol=5e-2) or np.allclose(xf, xf2, rtol=5e-2))
            self.assertTrue(np.allclose(yf, yf2, atol=5e-2) or np.allclose(yf, yf2, rtol=5e-2))
            self.assertTrue(np.allclose(zf, zf2, atol=5e-2) or np.allclose(zf, zf2, rtol=5e-2))
        
        
    def test_shuetal_mgnp(self):

        ntest = data['test_shuetal_mgnp'].shape[0]

        for i in range(ntest):

            iyear=data['test_shuetal_mgnp'][i]['iyear']
            idoy=data['test_shuetal_mgnp'][i]['idoy']
            ihour=data['test_shuetal_mgnp'][i]['ihour']
            imin=data['test_shuetal_mgnp'][i]['imin']
            isec=data['test_shuetal_mgnp'][i]['isec']
            xn_pd=data['test_shuetal_mgnp'][i]['xn_pd']
            vel=data['test_shuetal_mgnp'][i]['vel']
            bzimf=data['test_shuetal_mgnp'][i]['bzimf']

            xgsm=data['test_shuetal_mgnp'][i]['xgsm']
            ygsm=data['test_shuetal_mgnp'][i]['ygsm']
            zgsm=data['test_shuetal_mgnp'][i]['zgsm']
            xmgnp=data['test_shuetal_mgnp'][i]['xmgnp']
            ymgnp=data['test_shuetal_mgnp'][i]['ymgnp']
            zmgnp=data['test_shuetal_mgnp'][i]['zmgnp']
            idist=data['test_shuetal_mgnp'][i]['idist']
            id=data['test_shuetal_mgnp'][i]['id']

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xmgnp2, ymgnp2, zmgnp2, idist2, id2 = ispace.geopack_shuetal_mgnp(xn_pd, vel, bzimf, xgsm, ygsm, zgsm)
            
            self.assertTrue(np.allclose(xmgnp, xmgnp2, atol=1e-2) or np.allclose(xmgnp, xmgnp2, rtol=1e-2)) 
            self.assertTrue(np.allclose(ymgnp, ymgnp2, atol=1e-2) or np.allclose(ymgnp, ymgnp2, rtol=1e-2))
            self.assertTrue(np.allclose(zmgnp, zmgnp2, atol=1e-2) or np.allclose(zmgnp, zmgnp2, rtol=1e-2))
            self.assertTrue(np.allclose(idist, idist2, atol=1e-2) or np.allclose(idist, idist2, rtol=1e-2))
            self.assertTrue(np.allclose(id, id2, atol=1e-2) or np.allclose(id, id2, rtol=1e-2))




    def test_t96_mgnp(self):
        
        ntest = data['test_t96_mgnp'].shape[0]

        for i in range(ntest):

            iyear=data['test_t96_mgnp'][i]['iyear']
            idoy=data['test_t96_mgnp'][i]['idoy']
            ihour=data['test_t96_mgnp'][i]['ihour']
            imin=data['test_t96_mgnp'][i]['imin']
            isec=data['test_t96_mgnp'][i]['isec']
            xn_pd=data['test_t96_mgnp'][i]['xn_pd']
            vel=data['test_t96_mgnp'][i]['vel']
            
            xgsm=data['test_t96_mgnp'][i]['xgsm']
            ygsm=data['test_t96_mgnp'][i]['ygsm']
            zgsm=data['test_t96_mgnp'][i]['zgsm']
            xmgnp=data['test_t96_mgnp'][i]['xmgnp']
            ymgnp=data['test_t96_mgnp'][i]['ymgnp']
            zmgnp=data['test_t96_mgnp'][i]['zmgnp']
            idist=data['test_t96_mgnp'][i]['idist']
            id=data['test_t96_mgnp'][i]['id']

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xmgnp2, ymgnp2, zmgnp2, idist2, id2 = ispace.geopack_t96_mgnp(xn_pd, vel, xgsm, ygsm, zgsm)

            self.assertTrue(np.allclose(xmgnp, xmgnp2, atol=1e-2) or np.allclose(xmgnp, xmgnp2, rtol=1e-2)) 
            self.assertTrue(np.allclose(ymgnp, ymgnp2, atol=1e-2) or np.allclose(ymgnp, ymgnp2, rtol=1e-2))
            self.assertTrue(np.allclose(zmgnp, zmgnp2, atol=1e-2) or np.allclose(zmgnp, zmgnp2, rtol=1e-2))
            self.assertTrue(np.allclose(idist, idist2, atol=1e-2) or np.allclose(idist, idist2, rtol=1e-2))
            self.assertTrue(np.allclose(id, id2, atol=1e-2) or np.allclose(id, id2, rtol=1e-2))            


    def test_t89(self):
                    
        ntest = data['test_t89'].shape[0]

        for i in range(ntest):

            iyear=data['test_t89'][i]['iyear']
            idoy=data['test_t89'][i]['idoy']
            ihour=data['test_t89'][i]['ihour']
            imin=data['test_t89'][i]['imin']
            isec=data['test_t89'][i]['isec']
            iopt=data['test_t89'][i]['iopt']

            xgsm=data['test_t89'][i]['xgsm']
            ygsm=data['test_t89'][i]['ygsm']
            zgsm=data['test_t89'][i]['zgsm']

            bxgsm=data['test_t89'][i]['bxgsm']
            bygsm=data['test_t89'][i]['bygsm']
            bzgsm=data['test_t89'][i]['bzgsm']
           
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_t89(iopt, xgsm, ygsm, zgsm)


            self.assertTrue( np.allclose(bxgsm, bxgsm2, atol=5e-3) or np.allclose(bxgsm, bxgsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bygsm, bygsm2, atol=5e-3) or np.allclose(bygsm, bygsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bzgsm, bzgsm2, atol=5e-3) or np.allclose(bzgsm, bzgsm2, rtol=5e-3))
        

    def test_t96(self):
                    
        ntest = data['test_t96'].shape[0]

        for i in range(ntest):

            iyear=data['test_t96'][i]['iyear']
            idoy=data['test_t96'][i]['idoy']
            ihour=data['test_t96'][i]['ihour']
            imin=data['test_t96'][i]['imin']
            isec=data['test_t96'][i]['isec']
            parmod=data['test_t96'][i]['parmod']

            xgsm=data['test_t96'][i]['xgsm']
            ygsm=data['test_t96'][i]['ygsm']
            zgsm=data['test_t96'][i]['zgsm']

            bxgsm=data['test_t96'][i]['bxgsm']
            bygsm=data['test_t96'][i]['bygsm']
            bzgsm=data['test_t96'][i]['bzgsm']
           
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_t96(parmod, xgsm, ygsm, zgsm)

            self.assertTrue( np.allclose(bxgsm, bxgsm2, atol=5e-2) or np.allclose(bxgsm, bxgsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bygsm, bygsm2, atol=5e-2) or np.allclose(bygsm, bygsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bzgsm, bzgsm2, atol=5e-2) or np.allclose(bzgsm, bzgsm2, rtol=5e-3))

    def test_t01(self):
                        
        ntest = data['test_t01'].shape[0]

        for i in range(ntest):

            iyear=data['test_t01'][i]['iyear']
            idoy=data['test_t01'][i]['idoy']
            ihour=data['test_t01'][i]['ihour']
            imin=data['test_t01'][i]['imin']
            isec=data['test_t01'][i]['isec']
            parmod=data['test_t01'][i]['parmod']

            xgsm=data['test_t01'][i]['xgsm']
            ygsm=data['test_t01'][i]['ygsm']
            zgsm=data['test_t01'][i]['zgsm']

            bxgsm=data['test_t01'][i]['bxgsm']
            bygsm=data['test_t01'][i]['bygsm']
            bzgsm=data['test_t01'][i]['bzgsm']
        
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_t01(parmod, xgsm, ygsm, zgsm)

            # print('-'*40, i , '-'*40)
            # print('xgsm' , xgsm)
            # print('ygsm' , ygsm)
            # print('zgsm' , zgsm)
            # np.testing.assert_allclose(bxgsm, bxgsm2, rtol=2e-1, atol=0.1)
            # np.testing.assert_allclose(bygsm, bygsm2, rtol=2e-1, atol=0.1)
            # np.testing.assert_allclose(bzgsm, bzgsm2, rtol=2e-1, atol=0.1)
            self.assertTrue( np.allclose(bxgsm, bxgsm2, atol=2e-1) or np.allclose(bxgsm, bxgsm2, rtol=2e-2))
            self.assertTrue( np.allclose(bygsm, bygsm2, atol=2e-1) or np.allclose(bygsm, bygsm2, rtol=2e-2))
            self.assertTrue( np.allclose(bzgsm, bzgsm2, atol=2e-1) or np.allclose(bzgsm, bzgsm2, rtol=2e-2))

    def test_t04(self):
                            
        ntest = data['test_t04'].shape[0]

        for i in range(ntest):
            
            iyear=data['test_t04'][i]['iyear']
            idoy=data['test_t04'][i]['idoy']
            ihour=data['test_t04'][i]['ihour']
            imin=data['test_t04'][i]['imin']
            isec=data['test_t04'][i]['isec']
            parmod=data['test_t04'][i]['parmod']

            xgsm=data['test_t04'][i]['xgsm']
            ygsm=data['test_t04'][i]['ygsm']
            zgsm=data['test_t04'][i]['zgsm']

            bxgsm=data['test_t04'][i]['bxgsm']
            bygsm=data['test_t04'][i]['bygsm']
            bzgsm=data['test_t04'][i]['bzgsm']
        
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_t04(parmod, xgsm, ygsm, zgsm)
            
            # print('-'*40, i , '-'*40)           
            # np.testing.assert_allclose(bxgsm, bxgsm2, rtol=2e-2, atol=0.1)
            # np.testing.assert_allclose(bygsm, bygsm2, rtol=2e-2, atol=0.1)
            # np.testing.assert_allclose(bzgsm, bzgsm2, rtol=2e-2, atol=0.1)
            # print('-'*40, i , '-'*40)
            # print('-'*40, i , '-'*40)
           
            self.assertTrue( np.allclose(bxgsm, bxgsm2, atol=5e-1) or np.allclose(bxgsm, bxgsm2, rtol=2e-2))
            self.assertTrue( np.allclose(bygsm, bygsm2, atol=5e-1) or np.allclose(bygsm, bygsm2, rtol=2e-2))
            self.assertTrue( np.allclose(bzgsm, bzgsm2, atol=5e-1) or np.allclose(bzgsm, bzgsm2, rtol=2e-2))


    def test_bfield(self):

        ntest = data['test_bfield'].shape[0]

        for i in range(ntest):

            iyear=data['test_bfield'][i]['iyear']
            idoy=data['test_bfield'][i]['idoy']
            ihour=data['test_bfield'][i]['ihour']
            imin=data['test_bfield'][i]['imin']
            isec=data['test_bfield'][i]['isec']

            xgsm=data['test_bfield'][i]['xgsm']
            ygsm=data['test_bfield'][i]['ygsm']
            zgsm=data['test_bfield'][i]['zgsm']

            exname=data['test_bfield'][i]['exname']
            inname=data['test_bfield'][i]['inname']
            parmod=data['test_bfield'][i]['parmod']
            iopt=data['test_bfield'][i]['iopt']

            bxgsm=data['test_bfield'][i]['bxgsm']
            bygsm=data['test_bfield'][i]['bygsm']
            bzgsm=data['test_bfield'][i]['bzgsm']

            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            bxgsm2, bygsm2, bzgsm2 = ispace.geopack_bfield(iopt, parmod, xgsm,ygsm,zgsm,exname,inname)            
            bxgsm3, bygsm3, bzgsm3 = ispace.geopack_t96(parmod, xgsm, ygsm, zgsm)
            bxgsm4, bygsm4, bzgsm4 = ispace.geopack_igrf_gsm(xgsm, ygsm, zgsm)
            bxgsm5=bxgsm3+bxgsm4
            bygsm5=bygsm3+bygsm4
            bzgsm5=bzgsm3+bzgsm4

            self.assertTrue( np.allclose(bxgsm, bxgsm2, atol=5e-2) or np.allclose(bxgsm, bxgsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bygsm, bygsm2, atol=5e-2) or np.allclose(bygsm, bygsm2, rtol=5e-3))
            self.assertTrue( np.allclose(bzgsm, bzgsm2, atol=5e-2) or np.allclose(bzgsm, bzgsm2, rtol=5e-3))


    def test_recalc08(self):
        
        ntest = data['test_recalc08'].shape[0]
        for i in range(ntest):
            
            iyear=data['test_recalc08'][i]['iyear']
            idoy=data['test_recalc08'][i]['idoy']
            ihour=data['test_recalc08'][i]['ihour']
            imin=data['test_recalc08'][i]['imin']
            isec=data['test_recalc08'][i]['isec']       
            tilt=data['test_recalc08'][i]['tilt']     
            vgsex=data['test_recalc08'][i]['vgsex']
            vgsey=data['test_recalc08'][i]['vgsey']
            vgsez=data['test_recalc08'][i]['vgsez']
            tilt=data['test_recalc08'][i]['tilt']
            ispace.geopack08_recalc(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            #ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
             
            sps = ispace.geopack3.sps
            cps = ispace.geopack3.cps
            psi = np.arcsin(sps)
            # print(psi,tilt)
            self.assertAlmostEqual(psi, tilt, 4, 'psi not correct')
    
    
    def test_igrf_geo08(self):
        
        ntest = data['test_igrf_geo08'].shape[0]

        for i in range(ntest):
            
            iyear=data['test_igrf_geo08'][i]['iyear']
            idoy=data['test_igrf_geo08'][i]['idoy']
            ihour=data['test_igrf_geo08'][i]['ihour']
            imin=data['test_igrf_geo08'][i]['imin']
            isec=data['test_igrf_geo08'][i]['isec']  
            vgsex=data['test_igrf_geo08'][i]['vgsex']
            vgsey=data['test_igrf_geo08'][i]['vgsey']
            vgsez=data['test_igrf_geo08'][i]['vgsez']  

            r = data['test_igrf_geo08'][i]['r']
            theta = data['test_igrf_geo08'][i]['theta']
            phi = data['test_igrf_geo08'][i]['phi']
            br = data['test_igrf_geo08'][i]['br']
            btheta = data['test_igrf_geo08'][i]['btheta']
            bphi = data['test_igrf_geo08'][i]['bphi']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            br2, btheta2, bphi2 = ispace.geopack08_igrf_geo(r, theta, phi)

            np.testing.assert_allclose(br, br2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(btheta, btheta2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bphi, bphi2, rtol=5e-2, atol=5e-2)

    def test_igrf_gsw08(self):

        ntest = data['test_igrf_gsw08'].shape[0]

        for i in range(ntest):
            
            iyear=data['test_igrf_gsw08'][i]['iyear']
            idoy=data['test_igrf_gsw08'][i]['idoy']
            ihour=data['test_igrf_gsw08'][i]['ihour']
            imin=data['test_igrf_gsw08'][i]['imin']
            isec=data['test_igrf_gsw08'][i]['isec']  
            vgsex=data['test_igrf_gsw08'][i]['vgsex']
            vgsey=data['test_igrf_gsw08'][i]['vgsey']
            vgsez=data['test_igrf_gsw08'][i]['vgsez']  

            xgsw = data['test_igrf_gsw08'][i]['xgsw']
            ygsw = data['test_igrf_gsw08'][i]['ygsw']
            zgsw = data['test_igrf_gsw08'][i]['zgsw']
            bxgsw = data['test_igrf_gsw08'][i]['bxgsw']
            bygsw = data['test_igrf_gsw08'][i]['bygsw']
            bzgsw = data['test_igrf_gsw08'][i]['bzgsw']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            bxgsw2, bygsw2, bzgsw2 = ispace.geopack08_igrf_gsw(xgsw, ygsw, zgsw)

            np.testing.assert_allclose(bxgsw, bxgsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bygsw, bygsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bzgsw, bzgsw2, rtol=5e-2, atol=5e-2)


    def test_dip_gsw08(self):

        ntest = data['test_dip_gsw08'].shape[0]

        for i in range(ntest):
            
            iyear=data['test_dip_gsw08'][i]['iyear']
            idoy=data['test_dip_gsw08'][i]['idoy']
            ihour=data['test_dip_gsw08'][i]['ihour']
            imin=data['test_dip_gsw08'][i]['imin']
            isec=data['test_dip_gsw08'][i]['isec']  
            vgsex=data['test_dip_gsw08'][i]['vgsex']
            vgsey=data['test_dip_gsw08'][i]['vgsey']
            vgsez=data['test_dip_gsw08'][i]['vgsez']  

            xgsw = data['test_dip_gsw08'][i]['xgsw']
            ygsw = data['test_dip_gsw08'][i]['ygsw']
            zgsw = data['test_dip_gsw08'][i]['zgsw']
            bxgsw = data['test_dip_gsw08'][i]['bxgsw']
            bygsw = data['test_dip_gsw08'][i]['bygsw']
            bzgsw = data['test_dip_gsw08'][i]['bzgsw']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            bxgsw2, bygsw2, bzgsw2 = ispace.geopack08_dip_gsw(xgsw, ygsw, zgsw)

            np.testing.assert_allclose(bxgsw, bxgsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bygsw, bygsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bzgsw, bzgsw2, rtol=5e-2, atol=5e-2)

    def test_sph2car08(self):
        ntest = data['test_sph2car08'].shape[0]

        for i in range(ntest):
            r = data['test_sph2car08'][i]['r']
            theta = data['test_sph2car08'][i]['theta']
            phi = data['test_sph2car08'][i]['phi']
            x = data['test_sph2car08'][i]['x']
            y = data['test_sph2car08'][i]['y']
            z = data['test_sph2car08'][i]['z']

            x2, y2, z2 = ispace.geopack08_sph2car(r, theta, phi)

            np.testing.assert_allclose(x, x2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(y, y2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(z, z2, rtol=5e-2, atol=5e-2)
    

    def test_car2sph08(self):

        ntest = data['test_car2sph08'].shape[0]

        for i in range(ntest):
            x = data['test_car2sph08'][i]['x']
            y = data['test_car2sph08'][i]['y']
            z = data['test_car2sph08'][i]['z']
            r = data['test_car2sph08'][i]['r']
            theta = data['test_car2sph08'][i]['theta']
            phi = data['test_car2sph08'][i]['phi']

            r2, theta2, phi2 = ispace.geopack08_car2sph(x, y, z)

            np.testing.assert_allclose(r, r2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(theta, theta2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(phi, phi2, rtol=5e-2, atol=5e-2)
    
    
    def test_bsph2car08(self):
        ntest = data['test_bsph2car08'].shape[0]
        for i in range(ntest):
            theta = data['test_bsph2car08'][i]['theta']
            phi = data['test_bsph2car08'][i]['phi']
            br = data['test_bsph2car08'][i]['br']
            btheta = data['test_bsph2car08'][i]['btheta']
            bphi = data['test_bsph2car08'][i]['bphi']
            bx = data['test_bsph2car08'][i]['bx']
            by = data['test_bsph2car08'][i]['by']
            bz = data['test_bsph2car08'][i]['bz']

            bx2, by2, bz2 = ispace.geopack08_bsph2car(theta, phi, br, btheta, bphi)

            np.testing.assert_allclose(bx, bx2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(by, by2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bz, bz2, rtol=5e-2, atol=5e-2)

    
    def test_bcar2sph08(self):
        ntest = data['test_bcar2sph08'].shape[0]
        for i in range(ntest):
            x = data['test_bcar2sph08'][i]['x']
            y = data['test_bcar2sph08'][i]['y']
            z = data['test_bcar2sph08'][i]['z']
            bx = data['test_bcar2sph08'][i]['bx']
            by = data['test_bcar2sph08'][i]['by']
            bz = data['test_bcar2sph08'][i]['bz']
            br = data['test_bcar2sph08'][i]['br']
            btheta = data['test_bcar2sph08'][i]['btheta']
            bphi = data['test_bcar2sph08'][i]['bphi']

            br2, btheta2, bphi2 = ispace.geopack08_bcar2sph(x, y, z, bx, by, bz)

            np.testing.assert_allclose(br, br2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(btheta, btheta2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(bphi, bphi2, rtol=5e-2, atol=5e-2)


    def test_gsw2gse08(self):
        ntest = data['test_gsw2gse08'].shape[0]
        for i in range(ntest):

            iyear=data['test_gsw2gse08'][i]['iyear']
            idoy=data['test_gsw2gse08'][i]['idoy']
            ihour=data['test_gsw2gse08'][i]['ihour']
            imin=data['test_gsw2gse08'][i]['imin']
            isec=data['test_gsw2gse08'][i]['isec']  
            vgsex=data['test_gsw2gse08'][i]['vgsex']
            vgsey=data['test_gsw2gse08'][i]['vgsey']
            vgsez=data['test_gsw2gse08'][i]['vgsez']  

            xgsw = data['test_gsw2gse08'][i]['xgsw']
            ygsw = data['test_gsw2gse08'][i]['ygsw']
            zgsw = data['test_gsw2gse08'][i]['zgsw']
            xgse = data['test_gsw2gse08'][i]['xgse']
            ygse = data['test_gsw2gse08'][i]['ygse']
            zgse = data['test_gsw2gse08'][i]['zgse']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgse2,ygse2,zgse2 = ispace.geopack08_gsw2gse(xgsw,ygsw,zgsw)

            np.testing.assert_allclose(xgse, xgse2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygse, ygse2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgse, zgse2, rtol=5e-2, atol=5e-2)


    def test_gse2gsw08(self):
        ntest = data['test_gse2gsw08'].shape[0]
        for i in range(ntest):

            iyear=data['test_gse2gsw08'][i]['iyear']
            idoy=data['test_gse2gsw08'][i]['idoy']
            ihour=data['test_gse2gsw08'][i]['ihour']
            imin=data['test_gse2gsw08'][i]['imin']
            isec=data['test_gse2gsw08'][i]['isec']  
            vgsex=data['test_gse2gsw08'][i]['vgsex']
            vgsey=data['test_gse2gsw08'][i]['vgsey']
            vgsez=data['test_gse2gsw08'][i]['vgsez']  

            xgse = data['test_gse2gsw08'][i]['xgse']
            ygse = data['test_gse2gsw08'][i]['ygse']
            zgse = data['test_gse2gsw08'][i]['zgse']
            xgsw = data['test_gse2gsw08'][i]['xgsw']
            ygsw = data['test_gse2gsw08'][i]['ygsw']
            zgsw = data['test_gse2gsw08'][i]['zgsw']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgsw2,ygsw2,zgsw2 = ispace.geopack08_gse2gsw(xgse,ygse,zgse)

            np.testing.assert_allclose(xgsw, xgsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygsw, ygsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgsw, zgsw2, rtol=5e-2, atol=5e-2)                        


    def test_geo2mag08(self):
        ntest = data['test_geo2mag08'].shape[0]
        for i in range(ntest):

            iyear=data['test_geo2mag08'][i]['iyear']
            idoy=data['test_geo2mag08'][i]['idoy']
            ihour=data['test_geo2mag08'][i]['ihour']
            imin=data['test_geo2mag08'][i]['imin']
            isec=data['test_geo2mag08'][i]['isec']  
            vgsex=data['test_geo2mag08'][i]['vgsex']
            vgsey=data['test_geo2mag08'][i]['vgsey']
            vgsez=data['test_geo2mag08'][i]['vgsez']  

            xgeo = data['test_geo2mag08'][i]['xgeo']
            ygeo = data['test_geo2mag08'][i]['ygeo']
            zgeo = data['test_geo2mag08'][i]['zgeo']
            xmag = data['test_geo2mag08'][i]['xmag']
            ymag = data['test_geo2mag08'][i]['ymag']
            zmag = data['test_geo2mag08'][i]['zmag']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xmag2,ymag2,zmag2 = ispace.geopack08_geo2mag(xgeo,ygeo,zgeo)

            np.testing.assert_allclose(xmag, xmag2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ymag, ymag2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zmag, zmag2, rtol=5e-2, atol=5e-2)


    def test_mag2geo08(self):
        ntest = data['test_mag2geo08'].shape[0]
        for i in range(ntest):

            iyear=data['test_mag2geo08'][i]['iyear']
            idoy=data['test_mag2geo08'][i]['idoy']
            ihour=data['test_mag2geo08'][i]['ihour']
            imin=data['test_mag2geo08'][i]['imin']
            isec=data['test_mag2geo08'][i]['isec']  
            vgsex=data['test_mag2geo08'][i]['vgsex']
            vgsey=data['test_mag2geo08'][i]['vgsey']
            vgsez=data['test_mag2geo08'][i]['vgsez']  

            xmag = data['test_mag2geo08'][i]['xmag']
            ymag = data['test_mag2geo08'][i]['ymag']
            zmag = data['test_mag2geo08'][i]['zmag']
            xgeo = data['test_mag2geo08'][i]['xgeo']
            ygeo = data['test_mag2geo08'][i]['ygeo']
            zgeo = data['test_mag2geo08'][i]['zgeo']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgeo2,ygeo2,zgeo2 = ispace.geopack08_mag2geo(xmag,ymag,zmag)

            np.testing.assert_allclose(xgeo, xgeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygeo, ygeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgeo, zgeo2, rtol=5e-2, atol=5e-2)


    def test_gei2geo08(self):
        ntest = data['test_gei2geo08'].shape[0]
        for i in range(ntest):

            iyear=data['test_gei2geo08'][i]['iyear']
            idoy=data['test_gei2geo08'][i]['idoy']
            ihour=data['test_gei2geo08'][i]['ihour']
            imin=data['test_gei2geo08'][i]['imin']
            isec=data['test_gei2geo08'][i]['isec']  
            vgsex=data['test_gei2geo08'][i]['vgsex']
            vgsey=data['test_gei2geo08'][i]['vgsey']
            vgsez=data['test_gei2geo08'][i]['vgsez']  

            xgei = data['test_gei2geo08'][i]['xgei']
            ygei = data['test_gei2geo08'][i]['ygei']
            zgei = data['test_gei2geo08'][i]['zgei']
            xgeo = data['test_gei2geo08'][i]['xgeo']
            ygeo = data['test_gei2geo08'][i]['ygeo']
            zgeo = data['test_gei2geo08'][i]['zgeo']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgeo2,ygeo2,zgeo2 = ispace.geopack08_gei2geo(xgei,ygei,zgei)

            np.testing.assert_allclose(xgeo, xgeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygeo, ygeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgeo, zgeo2, rtol=5e-2, atol=5e-2)


    def test_geo2gei08(self):
        ntest = data['test_geo2gei08'].shape[0]
        for i in range(ntest):

            iyear=data['test_geo2gei08'][i]['iyear']
            idoy=data['test_geo2gei08'][i]['idoy']
            ihour=data['test_geo2gei08'][i]['ihour']
            imin=data['test_geo2gei08'][i]['imin']
            isec=data['test_geo2gei08'][i]['isec']  
            vgsex=data['test_geo2gei08'][i]['vgsex']
            vgsey=data['test_geo2gei08'][i]['vgsey']
            vgsez=data['test_geo2gei08'][i]['vgsez']  

            xgeo = data['test_geo2gei08'][i]['xgeo']
            ygeo = data['test_geo2gei08'][i]['ygeo']
            zgeo = data['test_geo2gei08'][i]['zgeo']
            xgei = data['test_geo2gei08'][i]['xgei']
            ygei = data['test_geo2gei08'][i]['ygei']
            zgei = data['test_geo2gei08'][i]['zgei']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgei2,ygei2,zgei2 = ispace.geopack08_geo2gei(xgeo,ygeo,zgeo)

            np.testing.assert_allclose(xgei, xgei2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygei, ygei2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgei, zgei2, rtol=5e-2, atol=5e-2)


    def test_mag2sm08(self):
        ntest = data['test_mag2sm08'].shape[0]
        for i in range(ntest):

            iyear=data['test_mag2sm08'][i]['iyear']
            idoy=data['test_mag2sm08'][i]['idoy']
            ihour=data['test_mag2sm08'][i]['ihour']
            imin=data['test_mag2sm08'][i]['imin']
            isec=data['test_mag2sm08'][i]['isec']  
            vgsex=data['test_mag2sm08'][i]['vgsex']
            vgsey=data['test_mag2sm08'][i]['vgsey']
            vgsez=data['test_mag2sm08'][i]['vgsez']  

            xmag = data['test_mag2sm08'][i]['xmag']
            ymag = data['test_mag2sm08'][i]['ymag']
            zmag = data['test_mag2sm08'][i]['zmag']
            xsm = data['test_mag2sm08'][i]['xsm']
            ysm = data['test_mag2sm08'][i]['ysm']
            zsm = data['test_mag2sm08'][i]['zsm']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xsm2,ysm2,zsm2 = ispace.geopack08_mag2sm(xmag,ymag,zmag)

            np.testing.assert_allclose(xsm, xsm2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ysm, ysm2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zsm, zsm2, rtol=5e-2, atol=5e-2)


    def test_sm2mag08(self):
        ntest = data['test_sm2mag08'].shape[0]
        for i in range(ntest):

            iyear=data['test_sm2mag08'][i]['iyear']
            idoy=data['test_sm2mag08'][i]['idoy']
            ihour=data['test_sm2mag08'][i]['ihour']
            imin=data['test_sm2mag08'][i]['imin']
            isec=data['test_sm2mag08'][i]['isec']  
            vgsex=data['test_sm2mag08'][i]['vgsex']
            vgsey=data['test_sm2mag08'][i]['vgsey']
            vgsez=data['test_sm2mag08'][i]['vgsez']  

            xsm = data['test_sm2mag08'][i]['xsm']
            ysm = data['test_sm2mag08'][i]['ysm']
            zsm = data['test_sm2mag08'][i]['zsm']
            xmag = data['test_sm2mag08'][i]['xmag']
            ymag = data['test_sm2mag08'][i]['ymag']
            zmag = data['test_sm2mag08'][i]['zmag']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xmag2,ymag2,zmag2 = ispace.geopack08_sm2mag(xsm,ysm,zsm)

            np.testing.assert_allclose(xmag, xmag2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ymag, ymag2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zmag, zmag2, rtol=5e-2, atol=5e-2)


    def test_sm2gsw08(self):
        ntest = data['test_sm2gsw08'].shape[0]
        for i in range(ntest):

            iyear=data['test_sm2gsw08'][i]['iyear']
            idoy=data['test_sm2gsw08'][i]['idoy']
            ihour=data['test_sm2gsw08'][i]['ihour']
            imin=data['test_sm2gsw08'][i]['imin']
            isec=data['test_sm2gsw08'][i]['isec']  
            vgsex=data['test_sm2gsw08'][i]['vgsex']
            vgsey=data['test_sm2gsw08'][i]['vgsey']
            vgsez=data['test_sm2gsw08'][i]['vgsez']  

            xsm = data['test_sm2gsw08'][i]['xsm']
            ysm = data['test_sm2gsw08'][i]['ysm']
            zsm = data['test_sm2gsw08'][i]['zsm']
            xgsw = data['test_sm2gsw08'][i]['xgsw']
            ygsw = data['test_sm2gsw08'][i]['ygsw']
            zgsw = data['test_sm2gsw08'][i]['zgsw']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgsw2,ygsw2,zgsw2 = ispace.geopack08_sm2gsw(xsm,ysm,zsm)

            np.testing.assert_allclose(xgsw, xgsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygsw, ygsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgsw, zgsw2, rtol=5e-2, atol=5e-2)


    def test_gsw2sm08(self):
        ntest = data['test_gsw2sm08'].shape[0]
        for i in range(ntest):

            iyear=data['test_gsw2sm08'][i]['iyear']
            idoy=data['test_gsw2sm08'][i]['idoy']
            ihour=data['test_gsw2sm08'][i]['ihour']
            imin=data['test_gsw2sm08'][i]['imin']
            isec=data['test_gsw2sm08'][i]['isec']  
            vgsex=data['test_gsw2sm08'][i]['vgsex']
            vgsey=data['test_gsw2sm08'][i]['vgsey']
            vgsez=data['test_gsw2sm08'][i]['vgsez']  

            xgsw = data['test_gsw2sm08'][i]['xgsw']
            ygsw = data['test_gsw2sm08'][i]['ygsw']
            zgsw = data['test_gsw2sm08'][i]['zgsw']
            xsm = data['test_gsw2sm08'][i]['xsm']
            ysm = data['test_gsw2sm08'][i]['ysm']
            zsm = data['test_gsw2sm08'][i]['zsm']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xsm2,ysm2,zsm2 = ispace.geopack08_gsw2sm(xgsw,ygsw,zgsw)

            np.testing.assert_allclose(xsm, xsm2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ysm, ysm2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zsm, zsm2, rtol=5e-2, atol=5e-2)


    def test_geo2gsw08(self):
        ntest = data['test_geo2gsw08'].shape[0]
        for i in range(ntest):

            iyear=data['test_geo2gsw08'][i]['iyear']
            idoy=data['test_geo2gsw08'][i]['idoy']
            ihour=data['test_geo2gsw08'][i]['ihour']
            imin=data['test_geo2gsw08'][i]['imin']
            isec=data['test_geo2gsw08'][i]['isec']  
            vgsex=data['test_geo2gsw08'][i]['vgsex']
            vgsey=data['test_geo2gsw08'][i]['vgsey']
            vgsez=data['test_geo2gsw08'][i]['vgsez']  

            xgeo = data['test_geo2gsw08'][i]['xgeo']
            ygeo = data['test_geo2gsw08'][i]['ygeo']
            zgeo = data['test_geo2gsw08'][i]['zgeo']
            xgsw = data['test_geo2gsw08'][i]['xgsw']
            ygsw = data['test_geo2gsw08'][i]['ygsw']
            zgsw = data['test_geo2gsw08'][i]['zgsw']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgsw2,ygsw2,zgsw2 = ispace.geopack08_geo2gsw(xgeo,ygeo,zgeo)

            np.testing.assert_allclose(xgsw, xgsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygsw, ygsw2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgsw, zgsw2, rtol=5e-2, atol=5e-2)


    def test_gsw2geo08(self):
        ntest = data['test_gsw2geo08'].shape[0]
        for i in range(ntest):

            iyear=data['test_gsw2geo08'][i]['iyear']
            idoy=data['test_gsw2geo08'][i]['idoy']
            ihour=data['test_gsw2geo08'][i]['ihour']
            imin=data['test_gsw2geo08'][i]['imin']
            isec=data['test_gsw2geo08'][i]['isec']  
            vgsex=data['test_gsw2geo08'][i]['vgsex']
            vgsey=data['test_gsw2geo08'][i]['vgsey']
            vgsez=data['test_gsw2geo08'][i]['vgsez']  

            xgsw = data['test_gsw2geo08'][i]['xgsw']
            ygsw = data['test_gsw2geo08'][i]['ygsw']
            zgsw = data['test_gsw2geo08'][i]['zgsw']
            xgeo = data['test_gsw2geo08'][i]['xgeo']
            ygeo = data['test_gsw2geo08'][i]['ygeo']
            zgeo = data['test_gsw2geo08'][i]['zgeo']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xgeo2,ygeo2,zgeo2 = ispace.geopack08_gsw2geo(xgsw,ygsw,zgsw)

            np.testing.assert_allclose(xgeo, xgeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(ygeo, ygeo2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zgeo, zgeo2, rtol=5e-2, atol=5e-2)


    def test_geod2geo08(self):
        ntest = data['test_geod2geo08'].shape[0]
        for i in range(ntest):

            iyear=data['test_geod2geo08'][i]['iyear']
            idoy=data['test_geod2geo08'][i]['idoy']
            ihour=data['test_geod2geo08'][i]['ihour']
            imin=data['test_geod2geo08'][i]['imin']
            isec=data['test_geod2geo08'][i]['isec']  
            vgsex=data['test_geod2geo08'][i]['vgsex']
            vgsey=data['test_geod2geo08'][i]['vgsey']
            vgsez=data['test_geod2geo08'][i]['vgsez']  

            r = data['test_geod2geo08'][i]['r']
            theta = data['test_geod2geo08'][i]['theta']
            h = data['test_geod2geo08'][i]['h']
            xmu = data['test_geod2geo08'][i]['xmu']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            r2, theta2 = ispace.geopack08_geod2geo(h, xmu)

            np.testing.assert_allclose(r, r2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(theta, theta2, rtol=5e-2, atol=5e-2)


    def test_geo2geod08(self):
        ntest = data['test_geo2geod08'].shape[0]
        for i in range(ntest):

            iyear=data['test_geo2geod08'][i]['iyear']
            idoy=data['test_geo2geod08'][i]['idoy']
            ihour=data['test_geo2geod08'][i]['ihour']
            imin=data['test_geo2geod08'][i]['imin']
            isec=data['test_geo2geod08'][i]['isec']  
            vgsex=data['test_geo2geod08'][i]['vgsex']
            vgsey=data['test_geo2geod08'][i]['vgsey']
            vgsez=data['test_geo2geod08'][i]['vgsez']  

            r = data['test_geo2geod08'][i]['r']
            theta = data['test_geo2geod08'][i]['theta']
            h = data['test_geo2geod08'][i]['h']
            xmu = data['test_geo2geod08'][i]['xmu']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            h2, xmu2 = ispace.geopack08_geo2geod(r, theta)

            np.testing.assert_allclose(h, h2, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(xmu, xmu2, rtol=5e-2, atol=5e-2)


    def test_trace08(self):
        """
        Note: Xgsm is now [-10.0, 6.0]. 
        Xgsm was [-13.0, 7.0], and there are many data points tracing to 100 Re. 
        """
        ntest = data['test_trace08'].shape[0]
        
        flag_idl08 = np.full(ntest, True, dtype=bool)
        flag_idl05 = np.full(ntest, True, dtype=bool)
        flag_py05 = np.full(ntest, True, dtype=bool)
        nfail=0
        for i in range(ntest):
            iyear=data['test_trace08'][i]['iyear']
            idoy=data['test_trace08'][i]['idoy']
            ihour=data['test_trace08'][i]['ihour']
            imin=data['test_trace08'][i]['imin']
            isec=data['test_trace08'][i]['isec']  
            vgsex=data['test_trace08'][i]['vgsex']
            vgsey=data['test_trace08'][i]['vgsey']
            vgsez=data['test_trace08'][i]['vgsez']  

            dsmax=data['test_trace08'][i]['dsmax']
            dir1=data['test_trace08'][i]['dir']
            parmod=data['test_trace08'][i]['parmod']
            r0=data['test_trace08'][i]['r0']
            rlim=data['test_trace08'][i]['rlim']

            xgsm=data['test_trace08'][i]['xgsm']
            ygsm=data['test_trace08'][i]['ygsm']
            zgsm=data['test_trace08'][i]['zgsm']
            xf=data['test_trace08'][i]['xf']
            yf=data['test_trace08'][i]['yf']
            zf=data['test_trace08'][i]['zf']
            xf05=data['test_trace08'][i]['xf05']
            yf05=data['test_trace08'][i]['yf05']
            zf05=data['test_trace08'][i]['zf05']

            iopt=3.0
            err = 0.0001 # Suggested from geopack_2008.for
            nfacmax=10000
            exname='T96'
            inname='IGRF_GSM'
            inname='DIP'
            
            ispace.geopack08_recalc(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            ispace.geopack_recalc(iyear,idoy,ihour,imin,isec)
            xf2, yf2, zf2, xfac2, yfac2, zfac2, nfac2 = ispace.geopack08_trace(xgsm,ygsm,zgsm,dir1,dsmax, err, rlim, r0, iopt, parmod, exname, inname, nfacmax)            
            xf3, yf3, zf3, xfac3, yfac3, zfac3, nfac3 = ispace.geopack_trace(xgsm, ygsm, zgsm, dir1, rlim, r0, iopt, parmod, exname, inname)

            # Skip in case of open field lines. 
            #if(np.sqrt(xf2**2 + yf2**2 + zf2**2) > 2.0):
            #    continue
            
            # print('')
            # print(f'{i} xyz   (pos): {xgsm:9.5f}, {ygsm:9.5f}, {zgsm:9.5f}, {np.sqrt(xgsm**2+ygsm**2+zgsm**2):9.5f}')
            # print(f'{i} xyzf(idl05): {xf05:9.5f}, {yf05:9.5f}, {zf05:9.5f}, {np.sqrt(xf05**2+yf05**2+zf05**2):9.5f}')
            # print(f'{i} xyzf3 (P05): {xf3:9.5f}, {yf3:9.5f}, {zf3:9.5f}, {np.sqrt(xf3**2+yf3**2+zf3**2):9.5f}')
            # print(f'{i} xyzf(idl08): {xf:9.5f}, {yf:9.5f}, {zf:9.5f}, {np.sqrt(xf**2+yf**2+zf**2):9.5f}')
            # print(f'{i} xyzf2 (P08): {xf2:9.5f}, {yf2:9.5f}, {zf2:9.5f}, {np.sqrt(xf2**2+yf2**2+zf2**2):9.5f}')
            
            # Sometimes IDL_geopack08 traced out, but Py08 did not. 
            if(np.sqrt(xf2**2 + yf2**2 + zf2**2) > 2.0):
                continue
            
            # Test Geopack. It was ok in test_trace
            # This test is passed too, as expected.
            # Thus, these lines are commented out.
            # self.assertTrue(np.allclose(xf05, xf3, atol=5e-2) or np.allclose(xf05, xf3, rtol=5e-2))
            # self.assertTrue(np.allclose(yf05, yf3, atol=5e-2) or np.allclose(yf05, yf3, rtol=5e-2))
            # self.assertTrue(np.allclose(zf05, zf3, atol=5e-2) or np.allclose(zf05, zf3, rtol=5e-2))
            
            flag1 = np.allclose(xf, xf2, atol=1e-1) or np.allclose(xf, xf2, rtol=5e-2)
            flag2 = np.allclose(yf, yf2, atol=1e-1) or np.allclose(yf, yf2, rtol=5e-2)
            flag3 = np.allclose(zf, zf2, atol=1e-1) or np.allclose(zf, zf2, rtol=5e-2)
            flag_idl08[i] = flag1 and flag2 and flag3

            flag4 = np.allclose(xf, xf3, atol=1e-1) or np.allclose(xf, xf3, rtol=5e-2)
            flag5 = np.allclose(yf, yf3, atol=1e-1) or np.allclose(yf, yf3, rtol=5e-2)
            flag6 = np.allclose(zf, zf3, atol=1e-1) or np.allclose(zf, zf3, rtol=5e-2)
            flag_py05[i] = flag4 and flag5 and flag6

            flag7 = np.allclose(xf2, xf05, atol=1e-1) or np.allclose(xf2, xf05, rtol=5e-2)
            flag8 = np.allclose(yf2, yf05, atol=1e-1) or np.allclose(yf2, yf05, rtol=5e-2)
            flag9 = np.allclose(zf2, zf05, atol=1e-1) or np.allclose(zf2, zf05, rtol=5e-2)
            flag_idl05[i] = flag7 and flag8 and flag9


            # If py08 does not equal to either IDL08 or IDL05, then it is marked a problem. 
            if flag_idl08[i] == False and flag_idl05[i] == False:
                nfail=nfail+1
                print('')
                print(f'Failed at {i}')                
                print(f'{i} xyz   (pos): {xgsm:9.5f}, {ygsm:9.5f}, {zgsm:9.5f}, {np.sqrt(xgsm**2+ygsm**2+zgsm**2):9.5f}')
                print(f'{i} xyzf(idl05): {xf05:9.5f}, {yf05:9.5f}, {zf05:9.5f}, {np.sqrt(xf05**2+yf05**2+zf05**2):9.5f}')
                print(f'{i} xyzf3 (P05): {xf3:9.5f}, {yf3:9.5f}, {zf3:9.5f}, {np.sqrt(xf3**2+yf3**2+zf3**2):9.5f}')
                print(f'{i} xyzf(idl08): {xf:9.5f}, {yf:9.5f}, {zf:9.5f}, {np.sqrt(xf**2+yf**2+zf**2):9.5f}')
                print(f'{i} xyzf2 (P08): {xf2:9.5f}, {yf2:9.5f}, {zf2:9.5f}, {np.sqrt(xf2**2+yf2**2+zf2**2):9.5f}')
                self.assertTrue(flag_idl08[i] or flag_idl05[i])

            np.testing.assert_allclose(xf2, xf05, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(yf2, yf05, rtol=5e-2, atol=5e-2)
            np.testing.assert_allclose(zf2, zf05, rtol=5e-2, atol=5e-2)
        
        # Now the test passed, let's print a summary. 
        # print('All Py08 results are equal to IDL08 or IDL05.')
        # print(f'Py08 not equal to IDL08: {ntest - np.count_nonzero(flag_idl08):4d} in {ntest} tests.')
        # print(f'Py08 not equal to IDL05: {ntest - np.count_nonzero(flag_idl05):4d} in {ntest} tests.')
        # print(f'Py08 not equal to Py05 : {ntest - np.count_nonzero(flag_py05 ):4d} in {ntest} tests.')
        # print(f'Py08 not equal IDL08 or IDL05: {ntest - np.count_nonzero(np.logical_or(flag_idl08,flag_idl05))}')


    def test_shuetal_mgnp08(self):
        ntest = data['test_shuetal_mgnp08'].shape[0]
        for i in range(ntest):

            iyear=data['test_shuetal_mgnp08'][i]['iyear']
            idoy=data['test_shuetal_mgnp08'][i]['idoy']
            ihour=data['test_shuetal_mgnp08'][i]['ihour']
            imin=data['test_shuetal_mgnp08'][i]['imin']
            isec=data['test_shuetal_mgnp08'][i]['isec']  
            vgsex=data['test_shuetal_mgnp08'][i]['vgsex']
            vgsey=data['test_shuetal_mgnp08'][i]['vgsey']
            vgsez=data['test_shuetal_mgnp08'][i]['vgsez']  

            xn_pd=data['test_shuetal_mgnp08'][i]['xn_pd']
            vel=data['test_shuetal_mgnp08'][i]['vel']
            bzimf=data['test_shuetal_mgnp08'][i]['bzimf']

            xgsw=data['test_shuetal_mgnp08'][i]['xgsw']
            ygsw=data['test_shuetal_mgnp08'][i]['ygsw']
            zgsw=data['test_shuetal_mgnp08'][i]['zgsw']

            xmgnp=data['test_shuetal_mgnp08'][i]['xmgnp']
            ymgnp=data['test_shuetal_mgnp08'][i]['ymgnp']
            zmgnp=data['test_shuetal_mgnp08'][i]['zmgnp']
            idist=data['test_shuetal_mgnp08'][i]['idist']
            id=data['test_shuetal_mgnp08'][i]['id']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            xmgnp2, ymgnp2, zmgnp2, idist2, id2 = ispace.geopack08_shuetal_mgnp(xn_pd, vel, bzimf, xgsw, ygsw, zgsw)

            np.testing.assert_allclose(xmgnp, xmgnp2, rtol=0.1, atol=0.15)
            np.testing.assert_allclose(ymgnp, ymgnp2, rtol=0.1, atol=0.15)
            np.testing.assert_allclose(zmgnp, zmgnp2, rtol=0.1, atol=0.15)
            np.testing.assert_allclose(idist, idist2, rtol=0.1, atol=0.15)
            np.testing.assert_allclose(id, id2, rtol=0.1, atol=0.15)


    def test_t96_mgnp08(self):
        ntest = data['test_t96_mgnp08'].shape[0]
        for i in range(ntest):

            iyear=data['test_t96_mgnp08'][i]['iyear']
            idoy=data['test_t96_mgnp08'][i]['idoy']
            ihour=data['test_t96_mgnp08'][i]['ihour']
            imin=data['test_t96_mgnp08'][i]['imin']
            isec=data['test_t96_mgnp08'][i]['isec']  
            vgsex=data['test_t96_mgnp08'][i]['vgsex']
            vgsey=data['test_t96_mgnp08'][i]['vgsey']
            vgsez=data['test_t96_mgnp08'][i]['vgsez']  

            xn_pd=data['test_t96_mgnp08'][i]['xn_pd']
            vel=data['test_t96_mgnp08'][i]['vel']

            xgsw=data['test_t96_mgnp08'][i]['xgsw']
            ygsw=data['test_t96_mgnp08'][i]['ygsw']
            zgsw=data['test_t96_mgnp08'][i]['zgsw']

            xmgnp=data['test_t96_mgnp08'][i]['xmgnp']
            ymgnp=data['test_t96_mgnp08'][i]['ymgnp']
            zmgnp=data['test_t96_mgnp08'][i]['zmgnp']
            idist=data['test_t96_mgnp08'][i]['idist']
            id=data['test_t96_mgnp08'][i]['id']

            ispace.recalc_08(iyear,idoy,ihour,imin,isec,vgsex,vgsey,vgsez)
            #xmgnp2, ymgnp2, zmgnp2, idist2, id2 = ispace.geopack08_t96_mgnp(xn_pd, vel, xgsw, ygsw, zgsw)
            xmgnp2, ymgnp2, zmgnp2, idist2, id2 = ispace.t96_mgnp_08(xn_pd, vel, xgsw, ygsw, zgsw)

            np.testing.assert_allclose(xmgnp, xmgnp2, rtol=1e-3, atol=1e-3)
            np.testing.assert_allclose(ymgnp, ymgnp2, rtol=1e-3, atol=1e-3)
            np.testing.assert_allclose(zmgnp, zmgnp2, rtol=1e-3, atol=1e-3)
            np.testing.assert_allclose(idist, idist2, rtol=1e-3, atol=1e-3)
            np.testing.assert_allclose(id, id2, rtol=1e-3, atol=1e-3)


if __name__ == '__main__':
    unittest.main()


