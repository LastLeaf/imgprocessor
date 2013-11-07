// Copyright (c) 2013 LastLeaf, MIT License, https://github.com/LastLeaf/imgprocessor
'use strict';

window.imgprocessor = (function(){

	// invert everything but alpha
	var revert = function(imgData) {
		var data = imgData.data;
		for(var i=0; i<data.length; i+=4) {
			data[i] = 255 - data[i];
			data[i+1] = 255 - data[i+1];
			data[i+2] = 255 - data[i+2];
		}
	};

	// seperate red
	var seperateRed = function(imgData) {
		var data = imgData.data;
		for(var i=0; i<data.length; i+=4) {
			data[i+1] = 0;
			data[i+2] = 0;
		}
	};

	// seperate green
	var seperateGreen = function(imgData) {
		var data = imgData.data;
		for(var i=0; i<data.length; i+=4) {
			data[i] = 0;
			data[i+2] = 0;
		}
	};

	// seperate blue
	var seperateBlue = function(imgData) {
		var data = imgData.data;
		for(var i=0; i<data.length; i+=4) {
			data[i] = 0;
			data[i+1] = 0;
		}
	};

	// monochrome
	var monochrome = function(imgData) {
		var data = imgData.data;
		for(var i=0; i<data.length; i+=4) {
			var t = 0.299*data[i] + 0.587*data[i+1] + 0.114*data[i+2];
			data[i]   = t;
			data[i+1] = t;
			data[i+2] = t;
		}
	}

	// greyness
	var setGrey = function(imgData, func){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var count = width*height;
		var y=0, u=0, v=0;
		for(var i=0; i<data.length; i+=4) {
			y = 0.299*data[i]+0.587*data[i+1]+0.114*data[i+2];
			u = 0.492*(data[i+2]-y[j]);
			v = 0.877*(data[i]-y[j]);
			y = func(y);
			data[i]   = y+1.140*v;
			data[i+1] = y-0.395*u-0.581*v;
			data[i+2] = y+2.032*u;
		}
	};

    // brightness
	var brightness = function(imgData, inc){
		setGrey(imgData, function(y){
			y += inc;
			if(y > 255) y = 255;
			if(y < 0) y = 0;
		});
	};

    // contrast
	var contrast = function(imgData, mul){
		setGrey(imgData, function(y){
			y = Math.round((y-128)*mul+128);
			if(y > 255) y = 255;
			if(y < 0) y = 0;
		});
	};

    // gamma correction
    var gammaCorrection = function(imgData, ratio){
		var data = imgData.data;
    	if(!ratio) ratio = 1;
    	var m = 10/ratio;
		for(var i=0; i<data.length; i+=4) {
			data[i] = Math.pow(data[i]/256, m) * 256;
			data[i+1] = Math.pow(data[i+1]/256, m) * 256;
			data[i+2] = Math.pow(data[i+2]/256, m) * 256;
		}
	};

    // histogram equalization
	var histogramEqualization = function(imgData) {
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var count = width*height;
		var y = new Float32Array(count);
		var u = new Float32Array(count);
		var v = new Float32Array(count);
		var percent = new Float32Array(256);
		var index = new uint32Array(256);
		var map = new uint32Array(256);
		for(var i=0,j=0;i<data.length;i+=4,j++) {
			y[j] = 0.299*data[i]+0.587*data[i+1]+0.114*data[i+2];
			u[j] = 0.492*(data[i+2]-y[j]);
			v[j] = 0.877*(data[i]-y[j]);				
		}
		for(var i=0; i<count; i++)
			index[Math.round(y[i])]++;
		percent[0] = index[0]/count*256;
		map[0] = Math.floor(percent[0]/2);
		for(var i=1; i<256; i++) {
			percent[i] = index[i]/count*256 + percent[i-1];
			map[i] = Math.floor((percent[i]+percent[i-1])/2);
		}
		for(var i=0;i<count;i++)
			y[i] = map[Math.round(y[i])];
		for(var i=0,j=0;i<data.length;i+=4,j++) {
			data[i]   = y[j]+1.140*v[j];
			data[i+1] = y[j]-0.395*u[j]-0.581*v[j];
			data[i+2] = y[j]+2.032*u[j];
		}
	};

	// FFT gauss low pass filter
	var lowPassGauss = function(fftData, d) {
		var real = fftData.real;
		var imag = fftData.imag;
		var width = fftData.width;
		var height = fftData.height;
		var xmid = width/2;
		var ymid = height/2;
		var t = -2*d*d;
		var p = 0;
		for(var y=0; y<height; y++) {
			for(var x=0; x<width; x++) {
				var k = Math.exp(((x-xmid)*(x-xmid)+(y-ymid)*(y-ymid))/t);
				real[p] *= k;
				imag[p] *= k;
				real[p+1] *= k;
				imag[p+1] *= k;
				real[p+2] *= k;
				imag[p+2] *= k;
				p += 4;
			}
		}
		return fftData;
	};

	// FFT gauss high pass filter
	var highPassGauss = function(fftData, d) {
		var real = fftData.real;
		var imag = fftData.imag;
		var width = fftData.width;
		var height = fftData.height;
		var xmid = width/2;
		var ymid = height/2;
		var t = -2*d*d;
		var p = 0;
		for(var y=0; y<height; y++) {
			for(var x=0; x<width; x++) {
				var k = 1 - Math.exp(((x-xmid)*(x-xmid)+(y-ymid)*(y-ymid))/t);
				real[p] *= k;
				imag[p] *= k;
				real[p+1] *= k;
				imag[p+1] *= k;
				real[p+2] *= k;
				imag[p+2] *= k;
				p += 4;
			}
		}
		return fftData;
	};

	// low pass filter
	var lowPass = function(imgData, arg) {
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var d = arg || 0;
		var fftData = CanvasFFT.FFT(data, width, height);
		CanvasFFT.IFFT(lowPassGauss(fftData, d), data, width, height);
	};

	// high pass filter
	var highPass = function(imgData, arg) {
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var d = arg || 0;
		var fftData = CanvasFFT.FFT(data, width, height);
		CanvasFFT.IFFT(highPassGauss(fftData, d), data, width, height);
	};

	// noice adder
	var noiceAdd = function(data, map) {
		var count = data.length/4;
		var Y = new Float32Array(count);
		var U = new Float32Array(count);
		var V = new Float32Array(count);
		for(var i=0, j=0; i<data.length; i+=4, j++) {
			Y[j] = 0.299*data[i]+0.587*data[i+1]+0.114*data[i+2];
			U[j] = 0.492*(data[i+2]-Y[j]);
			V[j] = 0.877*(data[i]-Y[j]);
		}
		for(var i=0; i<count; i++) {
			var t = Y[i] + map[ Math.floor(Math.random()*map.length) ];
			if(t<0) Y[i] = 0;
			else if(t>255) Y[i] = 255;
			else Y[i] = t;
		}
		for(var i=0, j=0; i<data.length; i+=4, j++) {
			data[i]   = Y[j]+1.140*V[j];
			data[i+1] = Y[j]-0.395*U[j]-0.581*V[j];
			data[i+2] = Y[j]+2.032*U[j];
		}
	};
	var noiceAddPd = function(data, pd, midValue) {
		var map = new Float32Array(256);
		var s = new Float32Array(256);
		s[0] = pd[0];
		for(var i=1; i<256; i++)
			s[i] = pd[i] + s[i-1];
		var m = 256/s[255];
		for(var i=0; i<256; i++)
			s[i] = Math.ceil(s[i]*m);
		var j = 0;
		for(var i=0; i<256; i++) {
			while(i >= s[j] && j < 255) j++;
			map[i] = j - midValue;
		}
		noiceAdd(data, map);
	};

	// gauss noice
	var noiceGauss = function(imgData, strength) {
		var data = imgData.data;
		var pd = new Float32Array(256);
		var r = Math.sqrt(strength);
		var PI = Math.PI;
		for(var i=0; i<256; i++)
			pd[i] = Math.exp(-i*i/(2*r*r)) / (Math.sqrt(2*PI)*r);
		noiceAddPd(data, pd, 0);
	};

	// uniform noice
	var noiceUniform = function(imgData, strength) {
		var data = imgData.data;
		var map = new Float32Array(255);
		var i=0;
		for(strength=Math.ceil(strength/2); strength; strength--) {
			map[i++] = -strength;
			map[i++] = strength;
		}
		for( ; i<256; i++)
			map[i] = 0;
		noiceAdd(data, map);
	};

	// two-value noice
	var noiceTwoValue = function(imgData, strength) {
		var data = imgData.data;
		var map = new Float32Array(256);
		var i=0;
		map[i++] = -strength;
		map[i++] = strength;
		for( ; i<256; i++)
			map[i] = 0;
		noiceAdd(data, map);
	};

	// mirror
	var mirror = function(imgData, isVertical){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		if(isVertical)
			for(var i=Math.floor(height/2)-1; i>=0; i--)
				for(var j=width-1; j>=0; j--) {
					var a = (i*width+j)*4;
					var b = ((height-1-i)*width+j)*4;
					var t = data[a];
					data[a] = data[b];
					data[b] = t;
					t = data[a+1];
					data[a+1] = data[b+1];
					data[b+1] = t;
					t = data[a+2];
					data[a+2] = data[b+2];
					data[b+2] = t;
					t = data[a+3];
					data[a+3] = data[b+3];
					data[b+3] = t;
				}
		else
			for(var i=height-1; i>=0; i--)
				for(var j=Math.floor(width/2)-1; j>=0; j--) {
					var a = (i*width+j)*4;
					var b = (i*width+(width-1-j))*4;
					var t = data[a];
					data[a] = data[b];
					data[b] = t;
					t = data[a+1];
					data[a+1] = data[b+1];
					data[b+1] = t;
					t = data[a+2];
					data[a+2] = data[b+2];
					data[b+2] = t;
					t = data[a+3];
					data[a+3] = data[b+3];
					data[b+3] = t;
				}
	};

	// emboss
	var emboss = function(imgData, depth){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		for(var i=height-1; i>=0; i--)
			for(var j=width-1; j>=0; j--) {
				var dest = (i*width+j)*4;
				if(i<depth || j<depth) {
					data[dest] = data[dest+1] = data[dest+2] = 128;
					continue;
				}
				var ori = ((i-depth)*width+j-depth)*4;
				data[dest]   -= data[ori] - 128;
				data[dest+1] -= data[ori+1] - 128;
				data[dest+2] -= data[ori+2] - 128;
			}
	};

	// generate guassian distribution
	var guassianDistribution = function(r2, rmax){
		var PI = Math.PI;
		// calc max ratio
		var r = 1;
		var min = 1/256 / (2*PI*r2);
		for(; r<=rmax; r++) {
			if(Math.exp(-r*r/(2*r2)) / (2*PI*r2) < min*(r*2-1)*(r*2-1))
				break;
		}
		r--;
		console.log(r);
		// generate guass
		var pds = r*2+1;
		var pd = new Float32Array(pds*pds);
		for(var i=-r; i<=r; i++)
			for(var j=-r; j<=r; j++) {
				var cur = (i+r)*pds + (j+r);
				pd[cur] = Math.exp(-(i*i+j*j)/(2*r2)) / (2*PI*r2);
			}
		var s = 0;
		for(var i=0; i<pd.length; i++)
			s += pd[i];
		for(var i=0; i<pd.length; i++)
			pd[i] /= s;
		return {
			r: r,
			size: pds,
			data: pd
		};
	};

	// gauss blur
	var blur = function(imgData, strength){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var gd = guassianDistribution(strength, (width>height?width:height));
		var rmax = gd.r;
		var size = gd.size;
		var pd = gd.data;
		// blur
		var blurColor = function(c){
			var t = new Uint32Array(width*height);
			for(var i=c, j=0; j<t.length; i+=4, j++)
				t[j] = data[i];
			for(var i=0; i<height; i++)
				for(var j=0; j<width; j++) {
					var dest = (i*width+j)*4+c;
					data[dest] = 0;
					for(var y=-rmax; y<=rmax; y++)
						for(var x=-rmax; x<=rmax; x++) {
							var cur = (y+rmax)*size + (x+rmax);
							if(i+y < 0 || i+y >= height || j+x < 0 || j+x >= width)
								var src = i*width+j;
							else
								var src = (i+y)*width + (j+x);
							data[dest] += pd[cur]*t[src];
						}
				}
		};
		blurColor(0);
		blurColor(1);
		blurColor(2);
	};

	// gauss shapen
	var shapen = function(imgData, strength){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var gd = guassianDistribution(strength, (width>height?width:height));
		var rmax = gd.r;
		var size = gd.size;
		var pd = gd.data;
		// blur
		var shapenColor = function(c){
			var t = new Uint32Array(width*height);
			for(var i=c, j=0; j<t.length; i+=4, j++)
				t[j] = data[i];
			for(var i=0; i<height; i++)
				for(var j=0; j<width; j++) {
					var dest = (i*width+j)*4+c;
					var r = data[dest]*2;
					for(var y=-rmax; y<=rmax; y++)
						for(var x=-rmax; x<=rmax; x++) {
							var cur = (y+rmax)*size + (x+rmax);
							if(i+y < 0 || i+y >= height || j+x < 0 || j+x >= width)
								var src = i*width+j;
							else
								var src = (i+y)*width + (j+x);
							r -= pd[cur]*t[src];
						}
					data[dest] = r;
				}
		};
		shapenColor(0);
		shapenColor(1);
		shapenColor(2);
	};

	// mozaic
	var mozaic = function(imgData, size){
		var data = imgData.data, width = imgData.width, height = imgData.height;
		var left = Math.floor(width%size/2);
		if(width < size) left = width;
		var top = Math.floor(height%size/2);
		if(height < size) top = height;
		// mean value function for area
		var mean = function(l, t, w, h){
			if(w <= 0 || h <= 0) return;
			var r=0, g=0, b=0, a=0;
			for(var i=0; i<h; i++)
				for(var j=0; j<w; j++) {
					var cur = ((i+t)*width + (j+l)) * 4;
					r += data[cur];
					g += data[cur+1];
					b += data[cur+2];
					a += data[cur+3];
				}
			r /= w*h;
			g /= w*h;
			b /= w*h;
			a /= w*h;
			for(var i=0; i<h; i++)
				for(var j=0; j<w; j++) {
					var cur = ((i+t)*width + (j+l)) * 4;
					data[cur] = r;
					data[cur+1] = g;
					data[cur+2] = b;
					data[cur+3] = a;
				}
		};
		// calc areas
		var calcRow = function(t, h){
			mean(0, t, left, h);
			for(var j=left; j<width-size; j+=size)
				mean(j, t, size, h);
			mean(j, t, width-j, h);
		};
		calcRow(0, top);
		for(var i=top; i<height-size; i+=size)
			calcRow(i, size);
		calcRow(i, height-i);
	};

	return {
		revert: revert,
		seperateRed: seperateRed,
		seperateGreen: seperateGreen,
		seperateBlue: seperateBlue,
		monochrome: monochrome,
		brightness: brightness,
		contrast: contrast,
		gammaCorrection: gammaCorrection,
		histogramEqualization: histogramEqualization,
		highPass: highPass,
		lowPass: lowPass,
		noiceGauss: noiceGauss,
		noiceUniform: noiceUniform,
		noiceTwoValue: noiceTwoValue,
		mirror: mirror,
		emboss: emboss,
		blur: blur,
		shapen: shapen,
		mozaic: mozaic
	};
})();

/* FFT for HTML5 canvas. Thanks the original and core author, Patrick Stein. http://nklein.com/2009/09/fourier-transforms-in-javascript/ */
var CanvasFFT = (function(){

	var __getNextPowerOfTwo = function( nn ) {
	    var pp = 1;
	    while ( pp < nn )
			pp *= 2;
	    return pp;
	};

	var __rearrangeSamples = function( _array, _offset, _ww, _stride ) {
	    var target = 0;
	    for ( var pos = 0; pos < _ww; ++pos ) {
		if ( target > pos ) {
		    for ( var kk = 0; kk < 4; ++kk ) {
			var tmp = _array[ target*_stride + kk + _offset ];
			_array[ target*_stride + kk + _offset ]
				= _array[ pos*_stride + kk + _offset ];
			_array[ pos*_stride + kk + _offset ] = tmp;
		    }
		}
		var mask = _ww;
		while ( ( target & ( mask >>= 1 ) ) ) {
		    target &= ~mask;
		}
		target |= mask;
	    }
	};

	var __shiftSamples = function( _samps, _base, _ww, _stride ) {
	    var mid = _base + _ww * _stride / 2;

	    for ( var ii=0; ii < _ww/2; ++ii ) {
		for ( var kk=0; kk < 3; ++kk ) {
		    var tmp = _samps[ _base + ii*_stride + kk ];
		    _samps[ _base + ii*_stride + kk ]
			= _samps[ mid + ii*_stride + kk ];
		    _samps[ mid + ii*_stride + kk ] = tmp;
		}
	    }
	};

	var __performFFT = function( _real, _imag, _ww, _hh, _dx, _dy, _inverse ) {
	    for ( var jj = 0; jj < _hh; ++jj ) {
			if ( _inverse ) {
			    __shiftSamples( _real, jj * _dy, _ww, _dx );
			    __shiftSamples( _imag, jj * _dy, _ww, _dx );
			}

			__rearrangeSamples( _real, jj * _dy, _ww, _dx );
			__rearrangeSamples( _imag, jj * _dy, _ww, _dx );

			var pi = 3.14159265358979323846264338327950288;
			var angularScale = ( _inverse ) ? pi : -pi;
			for ( var step = 1; step < _ww; step += step ) {
			    var delta = angularScale / step;
			    var sine = Math.sin( delta / 2.0 );
			    var fac_r = 1.0;
			    var fac_i = 0.0;
			    var mul_r = -2.0 * sine * sine;
			    var mul_i = Math.sin( delta );
			    for ( var group = 0; group < step; ++group ) {
					for ( var pair = group; pair < _ww; pair += step * 2 ) {
					    var match = pair + step;
					    for ( var kk = 0; kk < 3; ++kk ) {
							var p_index = jj * _dy + pair * _dx + kk;
							var m_index = jj * _dy + match * _dx + kk;
							var rr = _real[ m_index ];
							var ii = _imag[ m_index ];

							var prod_r = rr * fac_r - ii * fac_i;
							var prod_i = rr * fac_i + ii * fac_r;

							_real[ m_index ] = _real[ p_index ] - prod_r;
							_imag[ m_index ] = _imag[ p_index ] - prod_i;
							_real[ p_index ] += prod_r;
							_imag[ p_index ] += prod_i;
					    }
					}

					var inc_r = mul_r * fac_r - mul_i * fac_i;
					var inc_i = mul_r * fac_i + mul_i * fac_r;
					fac_r += inc_r;
					fac_i += inc_i;
			    }
			}

			if ( !_inverse ) {
			    __shiftSamples( _real, jj * _dy, _ww, _dx );
			    __shiftSamples( _imag, jj * _dy, _ww, _dx );
			}
	    }

	    return {
		    width: _ww,
		    height: _hh,
		    real: _real,
		    imag: _imag,
	    };
	};

	var createImageData = function(width, height) {
		var canvas = document.createElement('canvas');
		canvas.width = width;
		canvas.height = height;
		return canvas.getContext('2d').getImageData(0, 0, width, height);
	};

	var FFTImageData = function( fftData ) {
		var real = fftData.real;
		var imag = fftData.imag;

		var imageData = createImageData(fftData.width, fftData.height);
		var result = imageData.data;

	    for ( var pp=0; pp < result.length; pp += 4 ) {
		    var rr = real[ pp ];
		    var ii = imag[ pp ];
		    result[ pp ] = Math.sqrt( rr*rr + ii*ii );
		    var rr = real[ pp+1 ];
		    var ii = imag[ pp+1 ];
		    result[ pp+1 ] = Math.sqrt( rr*rr + ii*ii );
		    var rr = real[ pp+2 ];
		    var ii = imag[ pp+2 ];
		    result[ pp+2 ] = Math.sqrt( rr*rr + ii*ii );
			result[ pp+3 ] = 255;
	    }

	    return imageData;
	};

	var FFT = function( canvasData, width, height ) {
	    var ww = __getNextPowerOfTwo(width);
	    var hh = __getNextPowerOfTwo(height);

	    var real = new Float32Array(ww * hh * 4);
	    var imag = new Float32Array(ww * hh * 4);

	    var q = 0;
	    var p = 0;
	    for ( var j = 0; j<height; ++j ) {
	    	var i = 0;
		    for ( ; i<width*4; ++i ) {
				real[ q ] = canvasData[ p++ ] / 255.0;
				imag[ q++ ] = 0.0;
		    }
		    for ( ; i<ww*4; ++i ) {
				real[ q ] = 0.0;
				imag[ q++ ] = 0.0;
		    }
	    }
	    for( ; q<real.length; ++q) {
	    	real[q] = 0.0;
	    	imag[q] = 0.0;
	    }

	    var fftData = __performFFT( real, imag, ww, hh, 4, ww*4, false );
	    __performFFT( real, imag, hh, ww, ww*4, 4, false );

	    return fftData;
	};

	var IFFT = function( fftData, canvasData, width, height ) {
	    var ww = fftData.width;
	    var hh = fftData.height;
	    var real = fftData.real;
	    var imag = fftData.imag;

	    __performFFT( real, imag, hh, ww, ww*4, 4, true );
	    __performFFT( real, imag, ww, hh, 4, ww*4, true );

	    var scale = ww * hh;

	    var p = 0;
	    var q = 0;
	    for(var j=0; j<height; j++) {
	    	for(var i=0; i<width*4; i+=4) {
	    		for(var k=0; k<3; k++) {
	    			var t = real[q++] * 255.0 / scale;
	    			if(t<0) canvasData[p++] = 0;
	    			else if(t>255) canvasData[p++] = 255;
	    			else canvasData[p++] = t;
	    		}
	    		canvasData[p++] = 255;
	    		q++;
	    	}
	    	q += (ww-width)*4;
	    }

	    return {
	    	width: ww,
	    	height: hh,
	    	data: canvasData
	    };
	};

	return {
		FFT: FFT,
		IFFT: IFFT,
		FFTImageData: FFTImageData
	};
})();
