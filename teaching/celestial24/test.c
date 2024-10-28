{
	  NumberOfParticles = RotatingSphNumberOfParticles;
	  /*  --SIS Model-- 
	   *  M(r)   = Mfact*r                          ; where Mfact = 2sigma2/G
	   *  Pot(r) = Pfact*ln(r)                      ; where Pfact = 2sigma2
	   *  DF(E)  = DFact*[exp((Emax-eee)/sigma2)-1] ; where DFact = (2*M_PI*sigam2)^-1.5
	   *  Emax = Pot(rmax), eee = 0.5*v*v+pot(r) 
	   */
	  float Mfact, Pfact, DFact, sigma2, Emax;
	  float distf, Mtot, DMMeanDensity;
	  Mfact = SphereMass / SphereRadius;  /* Total mass in SIS DM assuming gas has the same distribution*/
	  Mtot = Mfact*0.5;                   /* 0.5 is maximum radius of the sphere in a simulation unit */
	  DMMeanDensity = (1.0-RotatingSphFgas)*Mtot/NumberOfParticles*(GridSize/GridVolume);
	  sigma2 = Mfact*GravitationalConstant/(8.0*M_PI);
	  Pfact  = 2.0*sigma2;
	  DFact  = 1.0/pow(2.0*M_PI*sigma2, 1.5);
	  Emax = Pfact*log(0.5);              /* pot(rmax)*/
	  printf("sigma2 for SIS = %f\n",sigma2);
	  /* Allocate space */ 
	  this->AllocateNewParticles(NumberOfParticles);
	  
	  //start distribute dm particles
	  //pts distribute
	  int npart, iter, itermax = 4000;
	  float r, vr, vt, vt1, vt2, pot;
	  float vmax, fmax, ppp, qqq, eee;
	  float phi, cost, sint, cosp, sinp, azi;
  
	  MT::MersenneTwist mtrnd; 
	  mtrnd.init_genrand(RotatingDM[2]); //initialize the Mersenne Twister. 
	  printf("start dm initialize distribute...\n");
	  for(npart = 0; npart < NumberOfParticles; npart++)
	    {
	      r = (Mtot/Mfact)*mtrnd.genrand_real1();
	      phi = 2.0*M_PI*mtrnd.genrand_real1();
	      cost = 2.0*(mtrnd.genrand_real1() - 0.5);
	      sint = sqrt(1-cost*cost);
	      cosp = cos(phi);
	      sinp = sin(phi);
	      
	      /* Compute the DM velocity */
          pot =  Pfact*log(r);
        
	      vmax = sqrt(2.0*fabs(Emax-pot));
	      fmax = DFact*exp((Emax-pot)/sigma2);
	      
	      iter = 0;
	      do{
		ppp = -2.0*cos(acos(mtrnd.genrand_real1())/3.0 - 2.0*M_PI/3.0);
		qqq = (1.0 - ppp*ppp)*mtrnd.genrand_real1();
		
		vr = vmax*ppp;
		vt = vmax*sqrt(qqq);
		eee = pot + 0.5*(vr*vr + vt*vt);
		distf = DFact*exp((Emax-eee)/sigma2);
		
		iter++;
	      } while((mtrnd.genrand_real1() > distf/fmax) && iter < itermax );
	      
	      if(iter >= itermax)
		printf("WARNING: verlocity interation over flow %d\n",iter);
	      
	      if (mtrnd.genrand_real1() >= 0.5) vr *= -1.0;
	      azi = 2.0*M_PI*mtrnd.genrand_real1();
	      vt1 = vt*cos(azi);
	      vt2 = vt*sin(azi);
		  		  
	      ParticleMass[npart] = DMMeanDensity;
	      ParticleNumber[npart] = npart;
	      ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
	      
	      ParticlePosition[0][npart] = 0.5 + r * sint*cosp;
	      ParticlePosition[1][npart] = 0.5 + r * sint*sinp;
	      ParticlePosition[2][npart] = 0.5 + r * cost;
	      
	      ParticleVelocity[0][npart] = vr * sint*cosp + vt1 * cost*cosp - vt2*sinp;
	      ParticleVelocity[1][npart] = vr * sint*sinp + vt1 * cost*sinp + vt2*cosp;
	      ParticleVelocity[2][npart] = vr * cost      - vt1 * sint;

		} /*End of NumberOfParticles*/
		printf("finishing dm distribution initialize ...\n");
	  
     }
