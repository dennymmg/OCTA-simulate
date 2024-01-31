// Angiography methods implemented
// Speckle Variance

clear
// helper functions to check the position of particles
exec('C:\Users\Denny\Documents\Research\angio\get_Particle_size.sci');

format("v",15)
scf(1)
clf

v1=0
v=(v1:50d-6:2d-3)

SIZE_P=(5*rand(1,6400)+3)*1d-6//randomly generating particles of the size between 3 to 8 um

i_ref=6d11
noise=3d4 // noise
i_samp=6d10
i_dc=i_ref+i_samp

ns=360 // no of pixels in camera
L1=840d-9
L2=860d-9
L=linspace(L1,L2,ns)
kk=2*%pi./L
k1=2*%pi/L1
k2=2*%pi/L2
k=linspace(k2,k1,ns)//extrapolating data in linear k.

fs=2*%pi/(k(2)-k(1))
df=(fs/ns)*0.5
depth1=50*df
z12=df  //interparticle separation distance
dt=10d-3 // 100 fps
M=40        // no. of frames 
texp=800d-6

SpeckVar=zeros(v)
nbp=depth1      // back position of first particle(We are assuming at the beggining that first particle is inside bin, where nbp is exactly touching depth1
for j=1:length(v)
    a=0;
    t=0;
    p=1                                            // particle index
    Fm=0
    Fmvec=zeros(M,1)
    phasediffvec=zeros(M,1)
    for i=1:M                                      //  M = no of frames 
        Size_p=SIZE_P(p)                           //original size of the particle at starting of each iteration(Size is randomy distributed between 6 to 10 micrometer
        Size=getParticleSize(depth1,df,nbp,Size_p) // size of particle at time 't'
        sk=i_dc+2*sqrt(i_ref*i_samp)*texp*cos(2*kk*(depth1+Size/2-v(j)*t-v(j)*texp/2)).*sinc(kk*v(j)*texp)*Size.*sinc(kk*Size)+noise*texp*rand(1,length(kk)) // Signal at time 't'
        nbp=nbp+v(j)*dt             //particle is moving v*dt amount of length in each iteration inside for loop.
        nfp=Size_p+nbp              //front position of of the particle
        t=t+dt
        s1=interp1(kk,sk,k,'linear') // converting data equispaced in lambda to data equispaced in k
        s1=s1//.*window('hn',length(s1))
        S1=fft(s1)
        s=S1(6:int(length(k)/2)) // igonring near zero frequencies and limiting f to fs/2
        [Fm1 Fm1i]=max(abs(s))
        Fmvec(i)=Fm1
        phase=atan(imag(s(Fm1i)),real(s(Fm1i)))
        if i~=1 // i not equal to 1
            a=a+2*Fm*Fm1/(Fm^2+Fm1^2)
        end
        if(isCompletelyRightOutside(depth1,df,nbp,Size))        //Current partcle has left the box completely. hence break and process new particle   
            nbp=nbp-z12+(SIZE_P(p)/2)-(SIZE_P((p+1)/2))         //set nbp of immediate next particle
            p=p+1                                               // changing particle index
        end
        Fm=Fm1
        Fmi=Fm1i // index of max value
        prevphase=phase
    end
    // Speckle Variance
    SpeckVar(j)=1/M*sum((Fmvec-1/M*sum(Fmvec)).^2)
end

plot(v,SpeckVar)
title('Speckle Variance','fontsize',5)
xlabel('Velocity --------->')
ylabel('Variance ----->')
