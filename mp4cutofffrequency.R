# MP4 max frequency clipping (FFmpeg encoding)
# www.overfitting.net
# https://www.overfitting.net/2024/03/ocultando-imagenes-en-el-espectro-de-un.html

library(tuneR)  # readWave() (read WAV audio files)

# Generate MP4 with white noise audio:
# ffmpeg -framerate 24 -i spectrogramgirl.png -i whitenoise_60s_48kHz.wav
# -c:v libx264 -crf 18 -pix_fmt yuv420p whitenoise.mp4

# Extract back audio data from MP4:
# ffmpeg -i whitenoise.mp4 whitenoise_60s_48kHz_afterMP4.wav


# Plot original vs MP4 white noise spectrum
plotres=400

sound=readWave("whitenoise_60s_48kHz.wav")  # read original WAV
waveform=sound@left
dft=abs(fft(waveform))
N=round(length(dft)/2)  # Primera mitad de la FFT
maxfreq=sound@samp.rate/2/1000  # Máx. frecuencia FFT en kHz
MAXDFT=max(dft)
values=dft[1:N]/MAXDFT

LEN=N/plotres
fftplot=array(0, plotres)
for (i in 1:plotres) {
    fftplot[i]=mean(values[((i-1)*LEN+1):(i*LEN)])
}
plot(seq(from=0, to=maxfreq, len=plotres), cex.axis=0.5,
     fftplot, main='FFT white noise original (blue) vs MP4 (red)', ylim=c(0,0.3),
     xlab='Frequency (kHz)', ylab='Amplitude (Lin.)', col='blue', type='l')
axis(side=1, at=c(0:maxfreq), cex.axis=0.5)
abline(h=mean(fftplot), col='lightgray')

sound=readWave("whitenoise_60s_48kHz_afterMP4.wav")  # read original WAV
waveform=sound@left
dft=abs(fft(waveform))
N=round(length(dft)/2)  # Primera mitad de la FFT
maxfreq=sound@samp.rate/2/1000  # Máx. frecuencia FFT en kHz
values=dft[1:N]/MAXDFT
LEN=N/plotres
fftplot=array(0, plotres)
for (i in 1:plotres) {
    fftplot[i]=mean(values[((i-1)*LEN+1):(i*LEN)])
}
lines(seq(from=0, to=maxfreq, len=plotres), fftplot, col='red', type='l')

