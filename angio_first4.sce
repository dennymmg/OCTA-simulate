// Angiography methods implemented
// 1. Intensity-based Differentiation
// 2. Intensity and phase based Differentiation 
// 3. Intensity-based Doppler Variance 
// 4. Phase-resolved Doppler Variance 

//helper functions to check the position of particles//
clear
exec('C:\Users\Denny\Documents\Research\angio\get_Particle_size.sci');

format("v",15)
scf(1)
clf

v1=0
v=(v1:50d-6:2d-3)

//randomly generating particles of the size between 3 to 8 um
SIZE_P=(5*rand(1,6400)+3)*1d-6

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
M=40     // no. of frames 
texp=800d-6

IbDiff=zeros(v)
IpbDiff=zeros(v)
IbDV=zeros(v)
PrDV=zeros(v)

nbp=depth1      // back position of first particle(We are assuming at the begining that first particle is inside bin, where nbp is exactly touching depth1
for j=1:length(v)
    a=0;b=0;c=0;d=0;e=0;f=0;
    t=0;
    p=1                                            // particle index
    Fm=0
    for i=1:M                                      //  M = no of frames 
        Size_p=SIZE_P(p)                           //original size of the particle at starting of each iteration
        Size=getParticleSize(depth1,df,nbp,Size_p) // size of particle at time 't'
        sk=i_dc+2*sqrt(i_ref*i_samp)*cos(2*kk*(depth1+Size/2-v(j)*t-v(j)*texp/2))*texp.*sinc(kk*v(j)*texp)*Size.*sinc(kk*Size)+noise*texp*rand(1,length(kk)) // Signal at time 't'
        nbp=nbp+v(j)*dt             //particle is moving v*dt amount of length in each iteration inside for loop.
        nfp=Size_p+nbp              //front position of of the particle
        t=t+dt
        s1=interp1(kk,sk,k,'linear') // converting data equispaced in lambda to data equispaced in k
        s1=s1//.*window('hn',length(s1))
        S1=fft(s1)
        s=S1(6:int(length(k)/2)) // igonring near zero frequencies and limiting f to fs/2
        [Fm1 Fm1i]=max(abs(s))
        Fm1=Fm1///texp
        if i~=1 // i not equal to 1
            a=a+Fm.*Fm1
            b=b+Fm^2
            c=c+Fm1^2
            d=d+abs(Fm1-Fm)// here s=Fm+1 and sprev=Fm
            e=e+sprev(Fmi)*conj(s(Fm1i))
            f=f+abs(sprev(Fmi)-s(Fm1i))
        end
        if(isCompletelyRightOutside(depth1,df,nbp,Size))        //Current partcle has left the box completely. hence break and process new particle   
            nbp=nbp-z12+(SIZE_P(p)/2)-(SIZE_P((p+1)/2))         //set nbp of immediate next particle
            p=p+1                                               // changing particle index
        end
        Fm=Fm1
        Fmi=Fm1i // index of max value
        sprev=s//texp
    end
    // Doppler variance Intensity-based
        IbDV(j)=(b+c-2*a)/(b+c)
    // Differentiation Intensity-based
    IbDiff(j)=(1/(M-1))*d
    // Doppler variance Phase-resolved
    PrDV(j)=(b+c-2*abs(e))/(b+c)
    // Differentiation Intensity and phase-based
    IpbDiff(j)=(1/(M-1))*f
end

subplot(2,2,1)
plot(v,IbDiff)
title('Intensity-based Differentiation','fontsize',5)
xlabel('Velocity --------->')
ylabel('Intensity ----->')
subplot(2,2,2)
plot(v,IpbDiff)
title('Intensity and phase based Differentiation','fontsize',5)
xlabel('Velocity --------->')
ylabel('Intensity ----->')
subplot(2,2,3)
plot(v,IbDV)
title('Intensity-based Doppler variance','fontsize',5)
xlabel('Velocity --------->')
ylabel('Variance ----->')
subplot(2,2,4)
plot(v,PrDV)
title('Phase-resolved Doppler variance','fontsize',5)
xlabel('Velocity --------->')
ylabel('Variance ----->')
