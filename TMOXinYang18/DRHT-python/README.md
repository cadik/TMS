This is the approximate implementation of the DRHT paper. The project page can be found here:
<p>https://ybsong00.github.io/cvpr18_imgcorrect/index.html</p>

The ldr2hdr model is from HDRCNN http://hdrv.org/hdrcnn/. We will release our original implementation of ldr2hdr soon. The model can be found at either the HDRCNN project page or here:
<p>https://drive.google.com/open?id=138JfKA5QzjDu78PLf6Ih9t5Qh2bBMvhI. </p> 
You need to download it at first and put it under checkpoint folder.

The illustration of the files and folders.

<p>############### folders ################</p>
<p>checkpoint      ---         pre-trained models</p>
<p>input           ---         input ldr images</p>
<p>hdr_output      ---         hdr files</p>
<p>samples         ---         ldr results</p>

<p>############### .py files ################</p>
<p>ldr2hdr.py and hdr2ldr.py define the ldr2hdr and hdr2ldr networks, respectively.</p>
<p>ldr2hdr_test.py and hdr2ldr_test.py provide simple evaluation.</p>

<p>############### notes ################</p>
<p>1. The ldr2hdr part is based on the Siggraph Asia 17 paper "HDR image reconstruction from a single exposure using deep CNNs".</p>
<p>2. The hdr2ldr part performs better when using large batch_size.</p>

<p>If you find the code useful, please cite the following papers:</p>

<pre><code>@inproceedings{yang-cvpr18-DRHT,
    author = {Yang, Xin and Xu, Ke and Song, Yibing and Zhang, Qiang and Wei, Xiaopeng and Rynson, Lau},
    title = {Image Correction via Deep Reciprocating HDR Transformation},
    booktitle = {IEEE Conference on Computer Vision and Pattern Recognition},
    year = {2018},
  }
</code></pre>

<pre><code>@article{EKDMU17,
  author       = "Eilertsen, Gabriel and Kronander, Joel, and Denes, Gyorgy and Mantiuk, Rafa\l and Unger, Jonas",
  title        = "HDR image reconstruction from a single exposure using deep CNNs",
  journal      = "ACM Transactions on Graphics (TOG)",
  year         = "2017",
}
</code></pre>
