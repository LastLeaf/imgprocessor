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

* `imgprocessor(_sourceImage_)` create an imgprocessor object. The _sourceImage_ can be an URL, <canvas>, <img>, or [ImageData](https://developer.mozilla.org/en-US/docs/Web/API/ImageData) object read from canvas.
* `.exec(_function(imageData){...}_)` run pending operations and call the callback function (if given). The callback function receives an ImageData object that can be put into canvases. Normally you should not modify the imageData.
* `.toCanvas(_function(img){...}_)` run pending operations and generate an <canvas>.
* `.toImage(_function(img){...}_)` run pending operations and generate an <img>. Currently, some browsers do not support it, so use `.toCanvas()` if possible.


Algorithm API
=============

The algorithm functions accepts at most 1 argument. Try the [demo](http://lastleaf.github.io/imgprocessor/) if you get confused with the following APIs.

* `.revert()`
* `.extractRed()`
* `.extractGreen()`
* `.extractBlue()`
* `.monochrome()`
* `.brightness(_addValue_)`
* `.contrast(_multiplyValue_)`
* `.gammaCorrection(_gamma_)`
* `.histogramEqualization()`
* `.noiceGauss(_strength_)`
* `.noiceUniform(_strength_)`
* `.noiceTwoValue(_strength_)`
* `.mirror(_isVertical_)`
* `.emboss(_depth_)`
* `.blur(_strength_)`
* `.shapen(_strength_)`
* `.mozaic(_blockSize_)`


LICENSE
=======

Copyright (c) 2013 LastLeaf, MIT License