Module Diagnostics_Custom
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Custom_MHD_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8  :: pfactor(my_r%min:my_r%max)
        Integer :: r,k, t
        pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                        /ref%density(my_r%min:my_r%max)
        
        
        
          !!!!!!!!!!!!!!!!!!!!!!!!! Advection Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        If (compute_quantity(curl_v_grad_v_r) .or. compute_quantity(curl_v_grad_v_r_squared)) Then
            If (compute_quantity(curl_v_grad_v_r)) Then
                DO_PSI
                    qty(PSI) = one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * DDBUFF(PSI,dvpdrdt) + &
                               buffer(PSI,dvrdt) * one_over_r(r) * buffer(PSI,dvpdr) + buffer(PSI,vtheta) * &
                               one_over_r(r) * DDBUFF(PSI,dvpdtdt) + buffer(PSI,dvtdt) * one_over_r(r) * & 
                               buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvpdtdp) - &
                               csctheta(t) * csctheta(t) * costheta(t) * buffer(PSI,dvpdp) + buffer(PSI,dvrdt) + &
                               cottheta(t) * buffer(PSI,dvtdt) - csctheta(t) * &
                               csctheta(t) * buffer(PSI,vtheta)) + buffer(PSI,dvpdt) * one_over_r(r) * (csctheta(t) * &
                               buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) + &
                               one_over_r(r) * costheta(t) * csctheta(t) * ref%density(r) * (buffer(PSI,vr) * &
                               one_over_r(r) * buffer(PSI,dvpdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                               buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * &
                               buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                               one_over_r(r) * csctheta(t) * ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvtdrdp) + &
                               buffer(PSI,dvrdp) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                               (DDBUFF(PSI,dvtdtdp) + buffer(PSI,dvrdp)) + buffer(PSI,dvtdp) * one_over_r(r) * &
                               (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * DDBUFF(PSI,dvtdpdp) - cottheta(t) * buffer(PSI,dvpdp)) + &
                               buffer(PSI,dvpdp) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvtdp) - &
                               cottheta(t) * buffer(PSI,vphi)))                                       
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! New squared + sqrt quantity (i.e. |curl_v_grad_v_r|)
            If (compute_quantity(curl_v_grad_v_r_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif            
            
        Endif
    
    
        If (compute_quantity(curl_v_grad_v_theta) .or. compute_quantity(curl_v_grad_v_theta_squared)) Then
            If (compute_quantity(curl_v_grad_v_theta)) Then
                DO_PSI
                    qty(PSI) = one_over_r(r) * csctheta(t) * ref%density(r) * ( buffer(PSI,vr) * DDBUFF(PSI,dvrdrdp) + &
                               buffer(PSI,dvrdp) * buffer(PSI,dvrdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                               (DDBUFF(PSI,dvrdtdp) - buffer(PSI,dvtdp)) + buffer(PSI,dvtdp) * one_over_r(r) * &
                               (buffer(PSI,dvrdt) - buffer(PSI,vtheta)) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * DDBUFF(PSI,dvrdpdp) - buffer(PSI,dvpdp)) + buffer(PSI,dvpdp) * &
                               one_over_r(r) * ( csctheta(t) * buffer(PSI,dvrdp) - buffer(PSI,vphi))) - &
                               one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * buffer(PSI,dvpdr) + &
                               buffer(PSI,vtheta) * one_over_r(r) * buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                               ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * DDBUFF(PSI,dvpdrdr) + buffer(PSI,dvrdr) * &
                               one_over_r(r) * buffer(PSI,dvpdr) - buffer(PSI,vr) * one_over_r(r) * one_over_r(r) * &
                               buffer(PSI,dvpdr) + buffer(PSI,vtheta) * one_over_r(r) * DDBUFF(PSI,dvpdrdt) + &
                               buffer(PSI,dvtdr) * one_over_r(r) * buffer(PSI,dvpdt) - buffer(PSI,vtheta) * one_over_r(r) * &
                               one_over_r(r) * buffer(PSI,dvpdt) + &
                               buffer(PSI,vtheta) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvpdrdp) + &
                               buffer(PSI,dvrdr) + cottheta(t) * buffer(PSI,dvtdr)) + &
                               buffer(PSI,dvpdr) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvpdp) + buffer(PSI,vr) + &
                               cottheta(t) * buffer(PSI,vtheta)) - buffer(PSI,vphi) *  one_over_r(r) * one_over_r(r) * &
                               (csctheta(t) *  buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                               ref%dlnrho(r) * &
                               ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * buffer(PSI,dvpdr) + &
                               buffer(PSI,vtheta) * one_over_r(r) * buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) 
                END_DO
                Call Add_Quantity(qty)
            Endif   
            If (compute_quantity(curl_v_grad_v_theta_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
        If (compute_quantity(curl_v_grad_v_phi) .or. compute_quantity(curl_v_grad_v_phi_squared)) Then
            If (compute_quantity(curl_v_grad_v_theta)) Then
                DO_PSI
                    qty(PSI) = ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvtdrdr) + buffer(PSI,dvrdr) * buffer(PSI,dvtdr) + &
                               buffer(PSI,vtheta) * one_over_r(r) * (DDBUFF(PSI,dvtdrdt) + buffer(PSI,dvrdr)) + &
                               buffer(PSI,dvtdr) * one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) - &
                               buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) + &
                               buffer(PSI,vtheta) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvtdrdp) - cottheta(t) * & 
                               buffer(PSI,dvpdr)) + buffer(PSI,dvpdr) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvtdp) - &
                               cottheta(t) * buffer(PSI,vphi)) - buffer(PSI,vphi) * one_over_r(r) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvtdp) - cottheta(t) * buffer(PSI,vphi))) + &   
                               ref%dlnrho(r) * ref%density(r) * (buffer(PSI,vr) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * &
                               one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvtdp) - cottheta(t) * buffer(PSI,vphi))) + &
                               one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * &
                               one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvtdp) - cottheta(t) * buffer(PSI,vphi))) - &
                               one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvrdrdt) + buffer(PSI,dvrdt) * &
                               buffer(PSI,dvrdr) + buffer(PSI,vtheta) * one_over_r(r) * (DDBUFF(PSI,dvrdtdt) + buffer(PSI,dvtdt)) + &
                               buffer(PSI,dvtdt) * one_over_r(r) * (buffer(PSI,dvrdt) - buffer(PSI,vtheta)) + &
                               buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvrdtdp) - csctheta(t) * csctheta(t) * &
                               costheta(t) * buffer(PSI,dvrdp) - buffer(PSI,dvpdt)) + buffer(PSI,dvpdt) * one_over_r(r) * &
                               (csctheta(t) * buffer(PSI,dvrdp) - buffer(PSI,vphi)))
                                  
                END_DO
                Call Add_Quantity(qty)
            Endif  
            If (compute_quantity(curl_v_grad_v_phi_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
                !!!!!!!!!!!!!!!!!!!!!!!!!! Buoyancy Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

        If (compute_quantity(curl_buoyancy_force_theta) .or. compute_quantity(curl_buoyancy_force_theta_squared)) Then
            If (compute_quantity(curl_buoyancy_force_theta)) Then
                DO_PSI
                    qty(PSI) = ref%Buoyancy_Coeff(r) * (csctheta(t) * &
                               radius(r) * buffer(PSI,dtdp))  ! since dtdp = (1/r)*d(temperature or entropy)/dphi
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_buoyancy_force_theta|)
            If (compute_quantity(curl_buoyancy_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        If (compute_quantity(curl_buoyancy_force_phi) .or. compute_quantity(curl_buoyancy_force_phi_squared)) Then
            If (compute_quantity(curl_buoyancy_force_phi)) Then
                DO_PSI
                    qty(PSI) = -ref%Buoyancy_Coeff(r) * ( radius(r) * buffer(PSI,dtdt))   ! since dtdt = (1/r)*d(temperature or $
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_buoyancy_force_phi|)
            If (compute_quantity(curl_buoyancy_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Magnetic Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_j_cross_b_r) .or. compute_quantity(curl_j_cross_b_r_squared)) Then
            If (compute_quantity(curl_j_cross_b_r)) Then
                DO_PSI
                    qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * buffer(PSI,dbtdt) * buffer(PSI,dbpdt) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdtdt) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * buffer(PSI,btheta) * &
                               buffer(PSI,bphi) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dbtdt) * &
                               buffer(PSI,bphi) + &                         
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * &
                               buffer(PSI,dbpdt) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * buffer(PSI,dbtdt) * buffer(PSI,dbtdp) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdtdp) + & 
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * &
                               buffer(PSI,dbtdp) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                               buffer(PSI,br) * buffer(PSI,dbrdp) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbrdt) * buffer(PSI,dbrdp) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdtdp) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,bphi) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdt) + &
                               one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,dbpdr) + &
                               one_over_r(r) * buffer(PSI,br) * DDBUFF(PSI,dbpdrdt) + &                             !(first part ended)
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               buffer(PSI,btheta) * buffer(PSI,dbpdt) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * costheta(t) * &
                               buffer(PSI,btheta) * buffer(PSI,bphi) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                               buffer(PSI,btheta) * buffer(PSI,dbtdp) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                               buffer(PSI,br) * buffer(PSI,dbrdp) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               buffer(PSI,br) * buffer(PSI,bphi) + &
                               one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,br) * buffer(PSI,dbpdr) - &    !(second part ended)
                               one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdp) - &
                               one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdp) - &
                               buffer(PSI,dbrdp) * buffer(PSI,dbtdr) - &
                               buffer(PSI,br) * DDBUFF(PSI,dbtdrdp) + &
                               buffer(PSI,dbrdp) * one_over_r(r) * buffer(PSI,dbrdt) + & 
                               buffer(PSI,br) * one_over_r(r) * DDBUFF(PSI,dbrdtdp) + &
                               one_over_r(r) * buffer(PSI,dbpdp) * buffer(PSI,dbpdt) + &
                               one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdtdp) + &
                               2 * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                               buffer(PSI,dbpdp) * buffer(PSI,dbtdp) - &               
                               buffer(PSI,bphi) * DDBUFF(PSI,dbtdpdp))
                END_DO
                Call Add_Quantity(qty)
            Endif
        
            ! New squared + sqrt quantity (i.e. |curl_j_cross_b_r|)
            If (compute_quantity(curl_j_cross_b_r_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif            
        
        Endif
        
        
        If (compute_quantity(curl_j_cross_b_theta) .or. compute_quantity(curl_j_cross_b_theta_squared)) Then
            If (compute_quantity(curl_j_cross_b_theta)) Then
                DO_PSI
                    qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                               buffer(PSI,dbpdp) * buffer(PSI,dbrdp) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                               buffer(PSI,bphi) * DDBUFF(PSI,dbrdpdp) - &          
                               2 * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                               one_over_r(r) * csctheta(t) * buffer(PSI,dbpdp) * buffer(PSI,dbpdr) - &
                               one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdp) - &
                               2 * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * buffer(PSI,dbtdp) - &    
                               one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbtdr) - &
                               one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbrdt) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdp) - &    !(first part ended)
                               one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               buffer(PSI,btheta) * buffer(PSI,bphi) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * buffer(PSI,dbtdp) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * &
                               buffer(PSI,br) * buffer(PSI,dbrdp) - &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                               one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) + &                            !(second part ended)
                               one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdrdt) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               buffer(PSI,btheta) * buffer(PSI,bphi) - &
                               one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dbtdr) * buffer(PSI,bphi) - &
                               one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * buffer(PSI,dbpdr) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI,dbtdr) * buffer(PSI,dbtdp) - &
                               buffer(PSI,btheta) * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbtdr) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * buffer(PSI,dbrdp) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI,dbrdr) * buffer(PSI,dbrdp) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdrdp) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                               one_over_r(r) * buffer(PSI,dbrdr) * buffer(PSI,bphi) - &
                               one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) - &
                               buffer(PSI,dbrdr) * buffer(PSI,dbpdr) - &
                               buffer(PSI,br) * DDBUFF(PSI,dbpdrdr))
                END_DO
                Call Add_Quantity(qty)
            Endif
        
            ! New squared + sqrt quantity (i.e. |curl_j_cross_b_theta|)
            If (compute_quantity(curl_j_cross_b_theta_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif
        
        
        If (compute_quantity(curl_j_cross_b_phi) .or. compute_quantity(curl_j_cross_b_phi_squared)) Then
            If (compute_quantity(curl_j_cross_b_phi)) Then
                DO_PSI
                    qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,br) + &
                               one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdr) - &
                               one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbrdt) - &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * (buffer(PSI, bphi))**2 + &
                               one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbtdp) - &     !(first part ended)
                               one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,br) + &
                               one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,br) + &
                               one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdr) + &
                               buffer(PSI,dbrdr) * buffer(PSI,dbtdr) + &
                               buffer(PSI,br) * DDBUFF(PSI,dbtdrdr) - &
                               buffer(PSI,dbrdr) * buffer(PSI,dbrdt) - &
                               one_over_r(r) * buffer(PSI,br) * DDBUFF(PSI,dbrdrdt) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,dbrdt) + &
                               one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * buffer(PSI,dbpdr) * buffer(PSI,dbpdt) - &
                               one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                               one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                               (buffer(PSI, bphi))**2 - &
                               2 * one_over_r(r) * csctheta(t) * costheta(t) * &
                               buffer(PSI, bphi) * buffer(PSI,dbpdr) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI, dbpdr) * buffer(PSI, dbtdp) - &
                               buffer(PSI, bphi) * one_over_r(r) * one_over_r(r) * csctheta(t)  * buffer(PSI, dbtdp) + &
                               one_over_r(r) * csctheta(t) * buffer(PSI, bphi) * DDBUFF(PSI,dbtdrdp) + &                 !(second part ended)
                               one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                               buffer(PSI, bphi) * buffer(PSI,dbrdp) - &
                               one_over_r(r) * csctheta(t) * buffer(PSI,dbpdt) * buffer(PSI,dbrdp) - &
                               one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbrdtdp) + & 
                               2 * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) + &
                               buffer(PSI,dbpdt) * buffer(PSI,dbpdr) + &
                               buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                               2 * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbtdt) + &
                               buffer(PSI,dbtdt) * buffer(PSI,dbtdr) + &
                               buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdt) - &
                               one_over_r(r) * buffer(PSI,dbtdt) * buffer(PSI,dbrdt) - &
                               one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdt))
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! New squared + sqrt quantity (i.e. |curl_j_cross_b_phi|)
            If (compute_quantity(curl_j_cross_b_phi_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!! Coriolis Force !!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_coriolis_force_r) .or. compute_quantity(curl_coriolis_force_r_squared)) Then
            If (compute_quantity(curl_coriolis_force_r)) Then
                DO_PSI
                    qty(PSI) = - ref%Coriolis_Coeff * ref%density(r) * one_over_r(r) * &
                               (-sintheta(t) * buffer(PSI,vtheta) + cottheta(t) * csctheta(t) * buffer(PSI,vtheta) + &
                               costheta(t) *  buffer(PSI,dvtdt) + 2 * costheta(t) *  buffer(PSI,vr) + sintheta(t) * &
                               buffer(PSI,dvrdt) + cottheta(t) *  buffer(PSI,dvpdt))    
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_coriolis_force_r|)
            If (compute_quantity(curl_coriolis_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif

        If (compute_quantity(curl_coriolis_force_theta) .or. compute_quantity(curl_coriolis_force_theta_squared)) Then
            If (compute_quantity(curl_coriolis_force_theta)) Then
                DO_PSI
                    qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (one_over_r(r) * (buffer(PSI,dvpdp) + &
                               costheta(t) *  buffer(PSI,vtheta) + sintheta(t) * buffer(PSI,vr)) + costheta(t) * &
                               buffer(PSI,dvtdr) + sintheta(t) * buffer(PSI,dvrdr) + ref%dlnrho(r) * costheta(t) * &
                               buffer(PSI,vtheta) + ref%dlnrho(r) * sintheta(t) * buffer(PSI,vr))
                END_DO
                Call Add_Quantity(qty)
            Endif
            
            ! New squared + sqrt quantity (i.e. |curl_coriolis_force_theta|)
            If (compute_quantity(curl_coriolis_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif
                
        Endif

        If (compute_quantity(curl_coriolis_force_phi) .or. compute_quantity(curl_coriolis_force_phi)) Then
            If (compute_quantity(curl_coriolis_force_phi)) Then
                DO_PSI
                    qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (ref%dlnrho(r) * costheta(t) * buffer(PSI,vphi) + &
                               costheta(t) * buffer(PSI,dvpdr) - one_over_r(r) * sintheta(t) * buffer(PSI,dvpdr)) 

                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_coriolis_force_phi|)
            If (compute_quantity(curl_coriolis_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif



        !!!!!!!!!!!!!!!!!!!!!!!!!!! Pressure Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       If (compute_quantity(curl_pressure_force_r) .or. compute_quantity(curl_pressure_force_r_squared)) Then
            If (compute_quantity(curl_pressure_force_r)) Then
                DO_PSI
                    qty(PSI) = pfactor(r) * &
                               OneOverRSquared(r) * csctheta(t) * (DDBUFF(PSI,dvpdtdp)-DDBUFF(PSI,dvpdtdp)) 
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_viscous_force_r|)
            If (compute_quantity(curl_pressure_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif           
        Endif
         
        If (compute_quantity(curl_pressure_force_theta) .or. compute_quantity(curl_pressure_force_theta_squared)) Then
            If (compute_quantity(curl_pressure_force_theta)) Then
                DO_PSI
                    qty(PSI) = pfactor(r) * &
                               OneOverRSquared(r) * csctheta(t) * (DDBUFF(PSI,dvpdrdp) + &
                               ref%dlnrho(r) * buffer(PSI,dpdp) + DDBUFF(PSI,dvpdrdp))
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_pressure_force_theta|)
            If (compute_quantity(curl_pressure_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif



        If (compute_quantity(curl_pressure_force_phi) .or. compute_quantity(curl_pressure_force_phi_squared)) Then
            If (compute_quantity(curl_pressure_force_phi)) Then
                DO_PSI
                    qty(PSI) = pfactor(r) * &
                               one_over_r(r) * csctheta(t) * (-DDBUFF(PSI,dvpdrdt) + DDBUFF(PSI,dvpdrdt) - &
                               ref%dlnrho(r) * buffer(PSI,dpdt))
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! New squared + sqrt quantity (i.e. |curl_viscous_force_phi|)
            If (compute_quantity(curl_pressure_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = sqrt(qty(PSI)*qty(PSI))
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif
        
 