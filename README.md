imgprocessor
============

An HTML5 image processing library. [The demo is here.](http://lastleaf.github.io/imgprocessor/)


How To Use It
=============

Include the following lines in your html files.

```html
<script type="text/javascript" src="src/imgprocessor_algorithm.js" imgprocessor-worker></script>
<script type="text/javascript" src="src/imgprocessor.js"></script>
```

The `imgprocessor-worker` attribute means that you need to try using [Web Worker](http://en.wikipedia.org/wiki/Web_worker). It allows processing in a seperate thread if browsers support.

Then you can process images in javascript code. The following code load an image from URL, blur it slightly, mirror it, and show it.

```js
imgprocessor('text.png').blur(2).mirror().toCanvas(function(canvas){
	document.body.appendChild(canvas);
});
```


Core API
========

* `imgprocessor( sourceImage )` create an imgprocessor object. The sourceImage can be an URL, <canvas>, <img>, or [ImageData](https://developer.mozilla.org/en-US/docs/Web/API/ImageData) object read from canvas.
* `.exec( function(imageData){...} )` run pending operations and call the callback function (if given). The callback function receives an ImageData object that can be put into canvases. Normally you should not modify the imageData.
* `.toCanvas( function(img){...} )` run pending operations and generate an <canvas>.
* `.toImage( function(img){...} )` run pending operations and generate an <img>. Currently, some browsers do not support it, so use `.toCanvas()` if possible.


Algorithm API
=============

The algorithm functions accepts at most 1 argument. Try the [demo](http://lastleaf.github.io/imgprocessor/) if you get confused with the following APIs.

* `.revert()`
* `.extractRed()`
* `.extractGreen()`
* `.extractBlue()`
* `.monochrome()`
* `.brightness( addValue )`
* `.contrast( multiplyValue )`
* `.gammaCorrection( gamma )`
* `.histogramEqualization()`
* `.noiceGauss( strength )`
* `.noiceUniform( strength )`
* `.noiceTwoValue( strength )`
* `.mirror( isVertical )`
* `.emboss( depth )`
* `.blur( strength )`
* `.shapen( strength )`
* `.mozaic( blockSize )`


LICENSE
=======

Copyright (c) 2013 LastLeaf, MIT License