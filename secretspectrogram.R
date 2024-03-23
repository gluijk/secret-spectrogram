# Spectrogram image encoding
# www.overfitting.net
# https://www.overfitting.net/2024/03/ocultando-imagenes-en-el-espectro-de-un.html

# Formas de codificar una imagen en el espectrograma de un sonido:
#  Generar una señal de ruido blanco y filtrarlo de forma adaptativa en el tiempo según la imagen a codificar
#  Considerar la imagen a codificar directamente como un espectrograma y ejecutar sobre él la FFT inversa a tramos temporales
#  Construir un generador multitono y modular sus distintas frecuencias en el tiempo en base a la imagen a codificar

library(tuneR)  # readWave() (read WAV audio file)
library(signal)  # specgram(), much faster than phonTools::spectrogram()
library(tiff)
library(png)
library(terra)  # resample()


# AUX FUNCTIONS

# Conversions time <-> samples. t=0 <-> n=1
time2sample=function(t, fs=48000) {(round(t*fs+1))}
sample2time=function(n, fs=48000) {((n-1)/fs)}

arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=DIMY, ncols=DIMX, extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method)
    return (as.array(rasterrs))
}


# SPECTROGRAM IMAGE ENCODING
fs=48000
fmin=500
fmax=20000
NFREQ=100  # number of total frequencies
DURATION=30  # duration in s

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
# matplot(t(signal), type='l')

image=readPNG("image.png")  # read image to encode
image=arrayresample(image, DIMX=ncol(signal), DIMY=nrow(signal))
image=matrix(image, nrow=nrow(image))
signalout=signal*image  # amplitude modulation of signal using image
rm(image, signal)  # free a lot of memory

signalout=colSums(signalout)  # add all frequency contributions
MAX=max(abs(signalout))
signalout=signalout/MAX

# Save .WAV
sound=readWave("sawtooth_bandlimited.wav")  # read arbitrary WAV

bits=16
sound@left=round(signalout*(2^(bits-1)-1))  # normalize to 16-bit integer
sound@samp.rate=fs
sound@bit=bits
writeWave(sound, filename="secretspectrogram.wav")



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
writeTIFF(P, "spectrogram.tif", bits.per.sample=16, compression="LZW")

# Photoshop processing (edited file saved as sRGB)
P=readTIFF("spectrogram2.tif")



# BUILD ANIMATION
datos=P
marco=readTIFF("marco.tif")  # frame of animation

DIMY=nrow(marco)
DIMX=ncol(marco)
CENTRO=DIMX/2
NFRAMES=ncol(datos)  # 4500frames, 30s audio track at 4500/30=150fps
WIDTH=80  # legend (title and frequencies)

# Create sequence
Offset=0
Offset2=0
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
    
    writePNG(frm, paste0("img", ifelse(frame<10, "000",
                                       ifelse(frame<100, "00",
                                              ifelse(frame<1000, "0", ""))), frame, ".png"))
    
    print(paste0(frame+1, "/", NFRAMES))
}


# Build MP4
# ffmpeg -framerate 150 -i img%4d.png -i secretspectrogram.wav -c:v libx264
#        -crf 18 -pix_fmt yuv420p secretspectrogram.mp4