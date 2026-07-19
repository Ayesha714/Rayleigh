Curl Momentum Equation
====================================================================

===================================================================================================================================================================== ====== ======================================= 
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r`                                                                               1301    curl\_v\_grad\_v\_r       
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta`                                                                          1302    curl\_v\_grad\_v\_theta   
 :math:`\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi`                                                                            1303    curl\_v\_grad\_v\_phi     
 :math:`\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r\right)^2`                                                                1304    curl\_v\_grad\_v\_r\_squared       
 :math:`\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta\right)^2`                                                           1305    curl\_v\_grad\_v\_theta\_squared   
 :math:`\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi\right)^2`                                                             1306    curl\_v\_grad\_v\_phi\_squared     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                                                                   1307    curl\_vp\_grad\_vm\_r     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`                                                              1308    curl\_vp\_grad\_vm\_theta 
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                                                                1309    curl\_vp\_grad\_vm\_phi   
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                                                                   1310    curl\_vm\_grad\_vp\_r     
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                                                              1311    curl\_vm\_grad\_vp\_theta 
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                                                1312    curl\_vm\_grad\_vp\_phi   
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r`                                                                             1313    curl\_vp\_grad\_vp\_r     
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta`                                                                        1314    curl\_vp\_grad\_vp\_theta 
 :math:`\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi`                                                                          1315    curl\_vp\_grad\_vp\_phi   
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r`                                                         1316    curl\_vm\_grad\_vm\_r     
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta`                                                    1317    curl\_vm\_grad\_vm\_theta 
 :math:`\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi`                                                      1318    curl\_vm\_grad\_vm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta`                                                                               1319    curl\_buoyancy\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi`                                                                                 1320    curl\_buoyancy\_force\_phi  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta\right)^2`                                                                1321    curl\_buoyancy\_force\_theta\_squared  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi\right)^2`                                                                  1322    curl\_buoyancy\_force\_phi\_squared  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\theta`                                                                    1323    curl\_buoyancy\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\phi`                                                                      1324    curl\_buoyancy\_pforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\theta`                                                                              1325    curl\_buoyancy\_mforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\phi`                                                                                1326    curl\_buoyancy\_mforce\_phi  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r`                                               1327    curl\_coriolis\_force\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta`                                          1328    curl\_coriolis\_force\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi`                                            1329    curl\_coriolis\_force\_phi    
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r\right)^2`                                1330    curl\_coriolis\_force\_r\_squared      
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta\right)^2`                           1331    curl\_coriolis\_force\_theta\_squared  
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi\right)^2`                             1332    curl\_coriolis\_force\_phi\_squared    
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_r`                                              1333    curl\_coriolis\_pforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\theta`                                         1334    curl\_coriolis\_pforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\phi`                                           1335    curl\_coriolis\_pforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_r`                              1336    curl\_coriolis\_mforce\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\theta`                         1337    curl\_coriolis\_mforce\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\phi`                           1338    curl\_coriolis\_mforce\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r`                                                   1339    curl\_viscous\_force\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta`                                              1340    curl\_viscous\_force\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi`                                                1341    curl\_viscous\_force\_phi     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_r`                                                  1342    curl\_viscous\_pforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\theta`                                             1343    curl\_viscous\_pforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\phi`                                               1344    curl\_viscous\_pforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_r`                                        1345    curl\_viscous\_mforce\_r      
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\theta`                                   1346    curl\_viscous\_mforce\_theta  
 :math:`\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\phi`                                     1347    curl\_viscous\_mforce\_phi    
 :math:`\left\nabla\times\nabla P\right]_r`                                                                                                                            1348    curl\_pressure\_force\_r       
 :math:`\left\nabla\times\nabla P\right]_\theta`                                                                                                                       1349    curl\_pressure\_force\_theta   
 :math:`\left\nabla\times\nabla P\right]_\phi`                                                                                                                         1350    curl\_pressure\_force\_phi     
 :math:`\left\nabla\times\nabla P\right]_r^2`                                                                                                                          1351    curl\_pressure\_force\_r\_squared       
 :math:`\left\nabla\times\nabla P\right]_\theta^2`                                                                                                                     1352    curl\_pressure\_force\_theta\_squared   
 :math:`\left\nabla\times\nabla P\right]_\phi^2`                                                                                                                       1353    curl\_pressure\_force\_phi\_squared     
 :math:`\left\nabla\times\nabla P\right]_r`                                                                                                                            1354    curl\_pressure\_pforce\_r      
 :math:`\left\nabla\times\nabla P\right]_\theta`                                                                                                                       1355    curl\_pressure\_pforce\_theta  
 :math:`\left\nabla\times\nabla P\right]_\phi`                                                                                                                         1356    curl\_pressure\_pforce\_phi    
 :math:`\left\nabla\times\nabla P\right]_r`                                                                                                                            1357    curl\_pressure\_mforce\_r      
 :math:`\left\nabla\times\nabla P\right]_\theta`                                                                                                                       1358    curl\_pressure\_mforce\_theta  
 :math:`\left\nabla\times\nabla P\right]_\phi`                                                                                                                         1359    curl\_pressure\_mforce\_phi    
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r`                              1360    curl\_j\_cross\_b\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta`                         1361    curl\_j\_cross\_b\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi`                           1362    curl\_j\_cross\_b\_phi   
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r\right)^2`               1363    curl\_j\_cross\_b\_r\_squared     
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta\right)^2`          1364    curl\_j\_cross\_b\_theta\_squared 
 :math:`\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi\right)^2`            1365    curl\_j\_cross\_b\_phi\_squared   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_r`                  1366    curl\_jp\_cross\_bm\_r       
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\theta`             1367    curl\_jp\_cross\_bm\_theta   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\phi`               1368    curl\_jp\_cross\_bm\_phi     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_r`                              1369    curl\_jm\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\theta`                         1370    curl\_jm\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\phi`                           1371    curl\_jm\_cross\_bp\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_r`                            1372    curl\_jm\_cross\_bm\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\theta`                       1373    curl\_jm\_cross\_bm\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\phi`                         1374    curl\_jm\_cross\_bm\_phi   
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_r`        1374    curl\_jp\_cross\_bp\_r     
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\theta`   1375    curl\_jp\_cross\_bp\_theta 
 :math:`\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\phi`     1376    curl\_jp\_cross\_bp\_phi   
===================================================================================================================================================================== ====== ======================================= 
