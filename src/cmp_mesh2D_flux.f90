! Copyright (c) 2020-2026 Damien Furfaro & Jacek Kosek
! SPDX-License-Identifier: LGPL-2.0-or-later

module cmp_mesh2D_flux_m
    use cmp_mesh2D_init_m
    implicit none

contains

subroutine required_variables_for_implicit_2D_heat_diffusion(me,nb_e,sim)
     
    type(mesh2D_t), intent(inout) :: me
    integer, intent(in) :: nb_e
    type(simulation_t), intent(in) :: sim

    real(dp) :: TL,TR,Tint,NormQTf,DerNormQTf,lambdL,dlambdLdT,lambdR,dlambdRdT,DerTintSELF,DerTintOTH,Fact,Num,Den
    integer :: jj, i
    logical :: foundEdgeLab

    TL=me%StVar%temp(nb_e)

    do i = 1, size(me%label%regionArr)
        if(trim(me%M2D_Prop%elem(nb_e)%PhysE)==trim(me%label%regionArr(i)%id)) then
            lambdL=me%label%regionArr(i)%mat%thermal_conductivity(TL)
            dlambdLdT=me%label%regionArr(i)%mat%thermal_conductivity_der(TL)
            exit
        endif
    enddo

    me%flxS%SomNormQT(nb_e)=0.0_dp
    me%flxS%DerSomNormQT(nb_e)=0.0_dp
    me%flxS%DerNormQTOth2D(nb_e,:)=0.0_dp
    do jj=1,3

        if((jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==2) .or. &
           (jj==3 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1) .or. &
           (jj==2 .and. me%M2D_Prop%elem(nb_e)%NbNeigh==1)) then       ! Boundary condition

            foundEdgeLab=.false.
            do i = 1, size(me%label%edgeArr)
                if(trim(me%M2D_Prop%elem(nb_e)%PhysN(jj))==trim(me%label%edgeArr(i)%id)) then
                    if(me%label%edgeArr(i)%type=='flux') then
                        NormQTf=-(me%label%edgeArr(i)%value%v0d()/(me%label%edgeArr(i)%LengthBC*me%extrusion_length))*&
                                me%M2D_Prop%elem(nb_e)%face(jj)
                        DerNormQTf=0.0_dp
                        foundEdgeLab=.true.
                        exit
                    else if(me%label%edgeArr(i)%type=='temp') then
                        NormQTf=-(lambdL/me%M2D_Prop%elem(nb_e)%Delta(jj))*(me%label%edgeArr(i)%value%v0d()-TL)*&
                                me%M2D_Prop%elem(nb_e)%face(jj)
                        DerNormQTf=-(me%M2D_Prop%elem(nb_e)%face(jj)/me%M2D_Prop%elem(nb_e)%Delta(jj))*&
                                   (dlambdLdT*(me%label%edgeArr(i)%value%v0d()-TL)-lambdL)       
                        foundEdgeLab=.true.
                        exit                 
                    else
                        print*,'edge label type not recognized - mesh2D'
                        stop
                    endif
                endif
            enddo
            if(.not.foundEdgeLab)then            ! SS Link expected
                NormQTf=0.0_dp
                DerNormQTf=0.0_dp
            endif

        else ! neighbor in the internal domain

            TR=me%StVar%temp(me%M2D_Prop%elem(nb_e)%Neigh(jj))

            do i = 1, size(me%label%regionArr)
                if(trim(me%M2D_Prop%elem(nb_e)%PhysN(jj))==trim(me%label%regionArr(i)%id)) then
                    lambdR=me%label%regionArr(i)%mat%thermal_conductivity(TR)
                    dlambdRdT=me%label%regionArr(i)%mat%thermal_conductivity_der(TR)
                    exit
                endif
            enddo

            Num=lambdL*me%M2D_Prop%elem(nb_e)%DeltaNeigh(jj)*TL+lambdR*me%M2D_Prop%elem(nb_e)%Delta(jj)*TR
            Den=lambdL*me%M2D_Prop%elem(nb_e)%DeltaNeigh(jj)+lambdR*me%M2D_Prop%elem(nb_e)%Delta(jj)
            Tint=Num/Den

            DerTintSELF=((dlambdLdT*me%M2D_Prop%elem(nb_e)%DeltaNeigh(jj)*TL+lambdL*me%M2D_Prop%elem(nb_e)%DeltaNeigh(jj))*Den-&
                        Num*(dlambdLdT*me%M2D_Prop%elem(nb_e)%DeltaNeigh(jj)))/(Den**2)   
            DerTintOTH=((dlambdRdT*me%M2D_Prop%elem(nb_e)%Delta(jj)*TR+lambdR*me%M2D_Prop%elem(nb_e)%Delta(jj))*Den-Num*&
                       (dlambdRdT*me%M2D_Prop%elem(nb_e)%Delta(jj)))/(Den**2)                  

            Fact=-me%M2D_Prop%elem(nb_e)%face(jj)/me%M2D_Prop%elem(nb_e)%Delta(jj)
            NormQTf=Fact*lambdL*(Tint-TL)

            DerNormQTf=Fact*(dlambdLdT*(Tint-TL)+lambdL*(DerTintSELF-1.0_dp))
            me%flxS%DerNormQTOth2D(nb_e,jj)=Fact*lambdL*DerTintOTH            

        endif

        me%flxS%SomNormQT(nb_e)=me%flxS%SomNormQT(nb_e)+NormQTf
        me%flxS%DerSomNormQT(nb_e)=me%flxS%DerSomNormQT(nb_e)+DerNormQTf
    enddo

end subroutine required_variables_for_implicit_2D_heat_diffusion   



subroutine heat_diffusion_mesh2D(me,sim)
    type(mesh2D_t), intent(inout) :: me
    type(simulation_t), intent(in) :: sim
    
    integer :: nb_e
    
    do nb_e=1,me%M2D_Prop%nb_elements
      call required_variables_for_implicit_2D_heat_diffusion(me,nb_e,sim)
    enddo
        
end subroutine heat_diffusion_mesh2D 

end module cmp_mesh2D_flux_m