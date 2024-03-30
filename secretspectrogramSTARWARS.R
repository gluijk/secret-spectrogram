# Spectrogram image encoding
# www.overfitting.net
# https://www.overfitting.net/2024/03/ocultando-imagenes-en-el-espectro-de-un.html

library(tuneR)  # readWave() (read WAV audio file)
library(signal)  # specgram(), much faster than phonTools::spectrogram()
library(tiff)
library(png)
library(terra)  # resample()


# AUX FUNCTIONS

# Conversions time <-> samples. t=0 <-> n=1
time2sample=function(t, fs=48000) {round(t*fs+1)}
sample2time=function(n, fs=48000) {(n-1)/fs}

# Input image resample
arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=DIMY, ncols=DIMX, extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method)
    return (as.array(rasterrs))
}



# SPECTROGRAM IMAGE ENCODING
fmin=17000  # min encoding frequency (Hz)
fmax=23000  # max encoding frequeny (Hz)
NFREQ=30  # total number of frequencies
DURATION=60+26.892  # generated audio clip duration (s)
fs=48000  # sampling frequency (Hz)
if (fmax>=fs/2) print("WARNING: to prevent aliasing fmax must be lower than fs/2")

type='lin'  # type=c('lin, 'log')
if (type=='lin') {
    freqs=seq(from=fmin, to=fmax, length.out=NFREQ)  # lin spaced freqs
} else {
    base=exp(1/(NFREQ-1)*log(fmax/fmin))
    freqs=fmin*base^(0:(NFREQ-1))  # log spaced freqs  
}

time=seq(from=0, to=DURATION, length.out=time2sample(DURATION))
signal=matrix(0, nrow=NFREQ, ncol=length(time))
for (i in 1:NFREQ) signal[NFREQ-i+1,]=sin(2*pi*freqs[i]*time)  # create tones
# Plot 10 lowest frequencies
samples1ms=time2sample(0.001)  # 1ms of signal
matplot(seq(from=0, to=1, length.out=samples1ms),
    t(signal[(NFREQ-9):NFREQ, 1:samples1ms]), type='l',
    main=paste0('10 lowest frequencies of ', NFREQ, ' frequencies'),
    xlab='Time (ms)', ylab='')
abline(h=0)


# LEFT CHANNEL
image=readPNG("imageSTARWARS.png")  # read image to encode. Font: "OPTINewsGothicLight"
image=arrayresample(image, DIMX=ncol(signal), DIMY=nrow(signal))
image=matrix(image, nrow=nrow(image))
signalout=signal*image  # amplitude modulation of signal using image
# rm(image, signal)  # free a lot of memory

signalout=colSums(signalout)  # add all frequency contributions
MAX=max(abs(signalout))
signalout=signalout/MAX


# RIGHT CHANNEL
image=readPNG("imageSTARWARS2.png")  # read image to encode. Font: "OPTINewsGothicLight"
image=arrayresample(image, DIMX=ncol(signal), DIMY=nrow(signal))
image=matrix(image, nrow=nrow(image))
signalout2=signal*image  # amplitude modulation of signal using image
rm(image, signal)  # free a lot of memory

signalout2=colSums(signalout2)  # add all frequency contributions
MAX=max(abs(signalout2))
signalout2=signalout2/MAX


# Save .WAV
sound=readWave("starwars.wav")  # read STAR WARS WAV

bits=16
signaladd=round(signalout*(2^(bits-1)-1))[1:length(sound@left)]  # normalize to 16-bit integer
sound@left=0.8*sound@left+0.2*signaladd

signaladd=round(signalout2*(2^(bits-1)-1))[1:length(sound@right)]  # normalize to 16-bit integer
sound@right=0.8*sound@right+0.2*signaladd

# sound@samp.rate=fs
# sound@bit=bits
sound=normalize(sound, unit="16")  # clip to normalize back to 16 bits
writeWave(sound, filename="secretspectrogramSTARWARS.wav")



# BUILD SPECTROGRAM
# Espectrograma
waveform=sound@left
waveform=waveform/max(abs(waveform))

# Create spectrogram
DIMY=512  # number of desired rows

# Spectrogram parameters
fftn=DIMY*2  # number of samples to use for the FFT
window=hanning(fftn)  # window size (samples)
overlap=ceiling(length(window)/2)  # standard overlap (samples)
overlap=ceiling(overlap*1.4)  # increase overlap to increase output ncol
step=fftn-overlap
print(paste0("FFT plot: nrow=", fftn/2, ", ncol=", round(length(waveform)/step)))
# nrow of FFT plot = fftn/2
# ncol of FFT plot = length(waveform)/step

# Calculate spectrogram
spec=specgram(x=waveform, n=fftn, Fs=fs, window=window, overlap=overlap)
P=abs(spec$S)  # discard phase information
P=P/max(P)  # normalize
P=P[nrow(P):1,]
hist(P, breaks=500)

# Log scaling and noise clipping
P=10*log(P)
NOISE=-90  # noise threshold in dB
P[P < NOISE]=NOISE
P=P-min(P)
P=P/max(P)

# Save spectrogram
writeTIFF(P, "spectrogramSTARWARS.tif", bits.per.sample=16, compression="LZW")

# Photoshop processing (edited file saved as sRGB)
P=readTIFF("spectrogram2STARWARS.tif")



# BUILD ANIMATION
datos=P
marco=readTIFF("marcoSTARWARS.tif")  # frame of animation with Darth Vader

DIMY=nrow(marco)
DIMX=ncol(marco)
CENTRO=DIMX/2
NFRAMES=ncol(datos)  # 21461frames, 86.892s audio track at 21461/86.892=246.98fps
WIDTH=80  # legend (title and frequencies)

# Create sequence
Offset=0
Offset2=0
frame2=0
for (frame in 0:(NFRAMES-1)) {
    frm=marco
    INI=CENTRO-frame
    FIN=DIMX
    
    if (INI<WIDTH+1) {
        Offset=Offset+1
        INI=WIDTH+1
    }
    INIREAD=1+Offset
    
    if (INIREAD+FIN-INI > NFRAMES) {
        Offset2=Offset2+1
        FIN=DIMX-Offset2
    }
    
    frm[,INI:FIN,]=datos[,INIREAD:(INIREAD+FIN-INI),]
    
    # Draw red line
    frm[,CENTRO,1]=1
    frm[,CENTRO,2:3]=0

    if (!frame%%8) {  # 2683frames, 86.892s audio track at 2683/86.89187=30.8775fps       
        writePNG(frm, paste0("img", ifelse(frame2<10, "0000",
                                    ifelse(frame2<100, "000",
                                    ifelse(frame2<1000, "00",
                                    ifelse(frame2<10000, "0", "")))),
                             frame2, ".png"))
        frame2=frame2+1
    }
    
    print(paste0(frame+1, "/", NFRAMES))
}



# BUILD MP4
# 247fps version:
# ffmpeg -framerate 246.98 -i img%5d.png -i secretspectrogramSTARWARS.wav -c:v libx264
#        -crf 18 -pix_fmt yuv420p secretspectrogramSTARWARS.mp4

# 30fps version:
# ffmpeg -framerate 30.8775 -i img%5d.png -i secretspectrogramSTARWARS.wav -c:v libx264
#        -crf 18 -pix_fmt yuv420p secretspectrogramSTARWARS.mp4