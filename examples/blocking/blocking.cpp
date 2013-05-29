  


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>fys4411/examples/Blocking/src/blocking.cpp at master · ComputationalPhysics/fys4411 · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="http://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <link rel="xhr-socket" href="/_sockets" />


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" />

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="YaNiN0Pb9SF6YbZceqiuybx2UFOv9AeI/iakB7nYRjQ=" name="csrf-token" />

    <link href="https://a248.e.akamai.net/assets.github.com/assets/github-aacfd01406222a8e32af6bf66a2eed1a08267178.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://a248.e.akamai.net/assets.github.com/assets/github2-30896c5685c3dd8da766a4fd3065a563107c9370.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://a248.e.akamai.net/assets.github.com/assets/frameworks-ec9348b8374c693b0749d0b95b215fe3f5414fd0.js" type="text/javascript"></script>
      <script src="https://a248.e.akamai.net/assets.github.com/assets/github-4844d00ec63aff537c875fbad136b46f55a18928.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="072f5e17a108af01465503ae7aeacf05">

        <link data-pjax-transient rel='permalink' href='/ComputationalPhysics/fys4411/blob/3f9d63eb6de6da4384f7d4b1ce1b0c2f5984dbbc/examples/Blocking/src/blocking.cpp'>
    <meta property="og:title" content="fys4411"/>
    <meta property="og:type" content="githubog:gitrepository"/>
    <meta property="og:url" content="https://github.com/ComputationalPhysics/fys4411"/>
    <meta property="og:image" content="https://secure.gravatar.com/avatar/5c6cbdbc94dfb05e77fd703de611d612?s=420&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png"/>
    <meta property="og:site_name" content="GitHub"/>
    <meta property="og:description" content="fys4411 - Repository for the course material for FYS4411."/>
    <meta property="twitter:card" content="summary"/>
    <meta property="twitter:site" content="@GitHub">
    <meta property="twitter:title" content="ComputationalPhysics/fys4411"/>

    <meta name="description" content="fys4411 - Repository for the course material for FYS4411." />


    <meta content="2341672" name="octolytics-dimension-user_id" /><meta content="7757310" name="octolytics-dimension-repository_id" />
  <link href="https://github.com/ComputationalPhysics/fys4411/commits/master.atom" rel="alternate" title="Recent Commits to fys4411:master" type="application/atom+xml" />

  </head>


  <body class="logged_out page-blob linux vis-public env-production  ">
    <div id="wrapper">

      
      
      

      
      <div class="header header-logged-out">
  <div class="container clearfix">

    <a class="header-logo-wordmark" href="https://github.com/">Github</a>

    <div class="header-actions">
      <a class="button primary" href="/signup">Sign up</a>
      <a class="button" href="/login?return_to=%2FComputationalPhysics%2Ffys4411%2Fblob%2Fmaster%2Fexamples%2FBlocking%2Fsrc%2Fblocking.cpp">Sign in</a>
    </div>

    <div class="command-bar js-command-bar  in-repository">


      <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
        <li class="features"><a href="/features">Features</a></li>
          <li class="enterprise"><a href="http://enterprise.github.com/">Enterprise</a></li>
          <li class="blog"><a href="/blog">Blog</a></li>
      </ul>
        <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">
  <a href="/search/advanced" class="advanced-search-icon tooltipped downwards command-bar-search" id="advanced_search" title="Advanced search"><span class="octicon octicon-gear "></span></a>

  <input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
      data-repo="ComputationalPhysics/fys4411"
      data-branch="master"
      data-sha="1f2a344ba4f73fe9e86ca5ea133bb4cf84d1ff52"
  >

    <input type="hidden" name="nwo" value="ComputationalPhysics/fys4411" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

  <div class="divider-vertical"></div>

</form>
    </div>

  </div>
</div>


      


            <div class="site hfeed" itemscope itemtype="http://schema.org/WebPage">
      <div class="hentry">
        
        <div class="pagehead repohead instapaper_ignore readability-menu ">
          <div class="container">
            <div class="title-actions-bar">
              

<ul class="pagehead-actions">



    <li>
      <a href="/login?return_to=%2FComputationalPhysics%2Ffys4411"
        class="minibutton js-toggler-target star-button entice tooltipped upwards"
        title="You must be signed in to use this feature" rel="nofollow">
        <span class="octicon octicon-star"></span>Star
      </a>
      <a class="social-count js-social-count" href="/ComputationalPhysics/fys4411/stargazers">
        5
      </a>
    </li>
    <li>
      <a href="/login?return_to=%2FComputationalPhysics%2Ffys4411"
        class="minibutton js-toggler-target fork-button entice tooltipped upwards"
        title="You must be signed in to fork a repository" rel="nofollow">
        <span class="octicon octicon-git-branch"></span>Fork
      </a>
      <a href="/ComputationalPhysics/fys4411/network" class="social-count">
        2
      </a>
    </li>
</ul>

              <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
                <span class="repo-label"><span>public</span></span>
                <span class="mega-octicon octicon-repo"></span>
                <span class="author vcard">
                  <a href="/ComputationalPhysics" class="url fn" itemprop="url" rel="author">
                  <span itemprop="title">ComputationalPhysics</span>
                  </a></span> /
                <strong><a href="/ComputationalPhysics/fys4411" class="js-current-repository">fys4411</a></strong>
              </h1>
            </div>

            
  <ul class="tabs">
    <li class="pulse-nav"><a href="/ComputationalPhysics/fys4411/pulse" class="js-selected-navigation-item " data-selected-links="pulse /ComputationalPhysics/fys4411/pulse" rel="nofollow"><span class="octicon octicon-pulse"></span></a></li>
    <li><a href="/ComputationalPhysics/fys4411" class="js-selected-navigation-item selected" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /ComputationalPhysics/fys4411">Code</a></li>
    <li><a href="/ComputationalPhysics/fys4411/network" class="js-selected-navigation-item " data-selected-links="repo_network /ComputationalPhysics/fys4411/network">Network</a></li>
    <li><a href="/ComputationalPhysics/fys4411/pulls" class="js-selected-navigation-item " data-selected-links="repo_pulls /ComputationalPhysics/fys4411/pulls">Pull Requests <span class='counter'>0</span></a></li>

      <li><a href="/ComputationalPhysics/fys4411/issues" class="js-selected-navigation-item " data-selected-links="repo_issues /ComputationalPhysics/fys4411/issues">Issues <span class='counter'>0</span></a></li>



    <li><a href="/ComputationalPhysics/fys4411/graphs" class="js-selected-navigation-item " data-selected-links="repo_graphs repo_contributors /ComputationalPhysics/fys4411/graphs">Graphs</a></li>


  </ul>
  
<div class="tabnav">

  <span class="tabnav-right">
    <ul class="tabnav-tabs">
          <li><a href="/ComputationalPhysics/fys4411/tags" class="js-selected-navigation-item tabnav-tab" data-selected-links="repo_tags /ComputationalPhysics/fys4411/tags">Tags <span class="counter blank">0</span></a></li>
    </ul>
  </span>

  <div class="tabnav-widget scope">


    <div class="select-menu js-menu-container js-select-menu js-branch-menu">
      <a class="minibutton select-menu-button js-menu-target" data-hotkey="w" data-ref="master">
        <span class="octicon octicon-branch"></span>
        <i>branch:</i>
        <span class="js-select-button">master</span>
      </a>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">

        <div class="select-menu-modal">
          <div class="select-menu-header">
            <span class="select-menu-title">Switch branches/tags</span>
            <span class="octicon octicon-remove-close js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-filters">
            <div class="select-menu-text-filter">
              <input type="text" id="commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
            </div>
            <div class="select-menu-tabs">
              <ul>
                <li class="select-menu-tab">
                  <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
                </li>
                <li class="select-menu-tab">
                  <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
                </li>
              </ul>
            </div><!-- /.select-menu-tabs -->
          </div><!-- /.select-menu-filters -->

          <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket css-truncate" data-tab-filter="branches">

            <div data-filterable-for="commitish-filter-field" data-filterable-type="substring">

                <div class="select-menu-item js-navigation-item selected">
                  <span class="select-menu-item-icon octicon octicon-check"></span>
                  <a href="/ComputationalPhysics/fys4411/blob/master/examples/Blocking/src/blocking.cpp" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="master" rel="nofollow" title="master">master</a>
                </div> <!-- /.select-menu-item -->
            </div>

              <div class="select-menu-no-results">Nothing to show</div>
          </div> <!-- /.select-menu-list -->


          <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket css-truncate" data-tab-filter="tags">
            <div data-filterable-for="commitish-filter-field" data-filterable-type="substring">

            </div>

            <div class="select-menu-no-results">Nothing to show</div>

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

  </div> <!-- /.scope -->

  <ul class="tabnav-tabs">
    <li><a href="/ComputationalPhysics/fys4411" class="selected js-selected-navigation-item tabnav-tab" data-selected-links="repo_source /ComputationalPhysics/fys4411">Files</a></li>
    <li><a href="/ComputationalPhysics/fys4411/commits/master" class="js-selected-navigation-item tabnav-tab" data-selected-links="repo_commits /ComputationalPhysics/fys4411/commits/master">Commits</a></li>
    <li><a href="/ComputationalPhysics/fys4411/branches" class="js-selected-navigation-item tabnav-tab" data-selected-links="repo_branches /ComputationalPhysics/fys4411/branches" rel="nofollow">Branches <span class="counter ">1</span></a></li>
  </ul>

</div>

  
  
  


            
          </div>
        </div><!-- /.repohead -->

        <div id="js-repo-pjax-container" class="container context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:05ee090275acf52daecc233f05ea3bda -->
<!-- blob contrib frag key: views10/v8/blob_contributors:v21:05ee090275acf52daecc233f05ea3bda -->


<div id="slider">
    <div class="frame-meta">

      <p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

        <div class="breadcrumb">
          <span class='bold'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ComputationalPhysics/fys4411" class="js-slide-to" data-branch="master" data-direction="back" itemscope="url"><span itemprop="title">fys4411</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ComputationalPhysics/fys4411/tree/master/examples" class="js-slide-to" data-branch="master" data-direction="back" itemscope="url"><span itemprop="title">examples</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ComputationalPhysics/fys4411/tree/master/examples/Blocking" class="js-slide-to" data-branch="master" data-direction="back" itemscope="url"><span itemprop="title">Blocking</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ComputationalPhysics/fys4411/tree/master/examples/Blocking/src" class="js-slide-to" data-branch="master" data-direction="back" itemscope="url"><span itemprop="title">src</span></a></span><span class="separator"> / </span><strong class="final-path">blocking.cpp</strong> <span class="js-zeroclipboard zeroclipboard-button" data-clipboard-text="examples/Blocking/src/blocking.cpp" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
        </div>

      <a href="/ComputationalPhysics/fys4411/find/master" class="js-slide-to" data-hotkey="t" style="display:none">Show File Finder</a>


        
  <div class="commit file-history-tease">
    <img class="main-avatar" height="24" src="https://secure.gravatar.com/avatar/556ca72a8c1d51dda3c9af642048bca8?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="24" />
    <span class="author"><a href="/sigvebs" rel="author">sigvebs</a></span>
    <time class="js-relative-date" datetime="2013-03-14T03:52:16-07:00" title="2013-03-14 03:52:16">March 14, 2013</time>
    <div class="commit-title">
        <a href="/ComputationalPhysics/fys4411/commit/3f9d63eb6de6da4384f7d4b1ce1b0c2f5984dbbc" class="message">Added blocking example.</a>
    </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong> contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2>Users on GitHub who have contributed to this file</h2>
      <ul class="facebox-user-list">
        <li>
          <img height="24" src="https://secure.gravatar.com/avatar/556ca72a8c1d51dda3c9af642048bca8?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="24" />
          <a href="/sigvebs">sigvebs</a>
        </li>
      </ul>
    </div>
  </div>


    </div><!-- ./.frame-meta -->

    <div class="frames">
      <div class="frame" data-permalink-url="/ComputationalPhysics/fys4411/blob/3f9d63eb6de6da4384f7d4b1ce1b0c2f5984dbbc/examples/Blocking/src/blocking.cpp" data-title="fys4411/examples/Blocking/src/blocking.cpp at master · ComputationalPhysics/fys4411 · GitHub" data-type="blob">

        <div id="files" class="bubble">
          <div class="file">
            <div class="meta">
              <div class="info">
                <span class="icon"><b class="octicon octicon-file-text"></b></span>
                <span class="mode" title="File Mode">file</span>
                  <span>100 lines (82 sloc)</span>
                <span>2.686 kb</span>
              </div>
              <div class="actions">
                <div class="button-group">
                      <a class="minibutton js-entice" href=""
                         data-entice="You must be signed in and on a branch to make or propose changes">Edit</a>
                  <a href="/ComputationalPhysics/fys4411/raw/master/examples/Blocking/src/blocking.cpp" class="button minibutton " id="raw-url">Raw</a>
                    <a href="/ComputationalPhysics/fys4411/blame/master/examples/Blocking/src/blocking.cpp" class="button minibutton ">Blame</a>
                  <a href="/ComputationalPhysics/fys4411/commits/master/examples/Blocking/src/blocking.cpp" class="button minibutton " rel="nofollow">History</a>
                </div><!-- /.button-group -->
              </div><!-- /.actions -->

            </div>
                <div class="blob-wrapper data type-c js-blob-data">
      <table class="file-code file-diff">
        <tr class="file-code-line">
          <td class="blob-line-nums">
            <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>

          </td>
          <td class="blob-line-code">
                  <div class="highlight"><pre><div class='line' id='LC1'><span class="cp">#include &quot;blocking.h&quot;</span></div><div class='line' id='LC2'><br/></div><div class='line' id='LC3'><span class="c1">//------------------------------------------------------------------------------</span></div><div class='line' id='LC4'><span class="n">Blocking</span><span class="o">::</span><span class="n">Blocking</span><span class="p">(</span><span class="n">string</span> <span class="n">dataPath</span><span class="p">)</span><span class="o">:</span></div><div class='line' id='LC5'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dataPath</span><span class="p">(</span><span class="n">dataPath</span><span class="p">)</span></div><div class='line' id='LC6'><span class="p">{</span></div><div class='line' id='LC7'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">deltaBlockSize</span> <span class="o">=</span> <span class="mi">10</span><span class="p">;</span></div><div class='line' id='LC8'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">nNodes</span> <span class="o">=</span> <span class="mi">4</span><span class="p">;</span></div><div class='line' id='LC9'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">maxBlockSizeTreshold</span> <span class="o">=</span> <span class="mf">1e4</span><span class="p">;</span></div><div class='line' id='LC10'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">stepLength</span> <span class="o">=</span> <span class="mi">10</span><span class="p">;</span></div><div class='line' id='LC11'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">outFilename</span> <span class="o">=</span> <span class="s">&quot;blocking.mat&quot;</span><span class="p">;</span></div><div class='line' id='LC12'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dataName</span> <span class="o">=</span> <span class="s">&quot;data&quot;</span><span class="p">;</span></div><div class='line' id='LC13'><span class="p">}</span></div><div class='line' id='LC14'><span class="c1">//------------------------------------------------------------------------------</span></div><div class='line' id='LC15'><span class="kt">void</span> <span class="n">Blocking</span><span class="o">::</span><span class="n">doBlocking</span><span class="p">()</span></div><div class='line' id='LC16'><span class="p">{</span></div><div class='line' id='LC17'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">cout</span> <span class="o">&lt;&lt;</span> <span class="s">&quot;Starting blocking analysis.&quot;</span> <span class="o">&lt;&lt;</span> <span class="n">endl</span><span class="p">;</span></div><div class='line' id='LC18'><br/></div><div class='line' id='LC19'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vec</span> <span class="n">data</span> <span class="o">=</span> <span class="n">readDataFromFile</span><span class="p">();</span></div><div class='line' id='LC20'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">N</span> <span class="o">=</span> <span class="n">data</span><span class="p">.</span><span class="n">n_elem</span><span class="p">;</span></div><div class='line' id='LC21'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">maxBlockSize</span> <span class="o">=</span> <span class="p">(</span><span class="kt">int</span><span class="p">)</span> <span class="n">N</span><span class="o">/</span><span class="mi">10</span><span class="p">;</span></div><div class='line' id='LC22'><br/></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">maxBlockSize</span> <span class="o">&gt;</span> <span class="n">maxBlockSizeTreshold</span><span class="p">)</span></div><div class='line' id='LC24'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">maxBlockSize</span> <span class="o">=</span> <span class="n">maxBlockSizeTreshold</span><span class="p">;</span></div><div class='line' id='LC25'><br/></div><div class='line' id='LC26'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">mat</span> <span class="nf">results</span><span class="p">(</span><span class="mi">3</span><span class="p">,</span> <span class="n">maxBlockSize</span><span class="o">/</span><span class="n">deltaBlockSize</span><span class="p">);</span></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">string</span> <span class="n">outFile</span> <span class="o">=</span> <span class="n">dataPath</span> <span class="o">+</span> <span class="s">&quot;/&quot;</span> <span class="o">+</span> <span class="n">outFilename</span><span class="p">;</span></div><div class='line' id='LC28'><br/></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">i</span><span class="o">*</span><span class="n">deltaBlockSize</span> <span class="o">&lt;=</span> <span class="n">maxBlockSize</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">blockSize</span> <span class="o">=</span> <span class="n">i</span><span class="o">*</span><span class="n">deltaBlockSize</span><span class="p">;</span></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">results</span><span class="p">.</span><span class="n">col</span><span class="p">(</span><span class="n">i</span><span class="o">-</span><span class="mi">1</span><span class="p">)</span> <span class="o">=</span> <span class="n">computeBlock</span><span class="p">(</span><span class="n">data</span><span class="p">,</span> <span class="n">blockSize</span><span class="p">);</span></div><div class='line' id='LC32'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC33'><br/></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">results</span><span class="p">.</span><span class="n">save</span><span class="p">(</span><span class="n">outFile</span><span class="p">,</span> <span class="n">arma_ascii</span><span class="p">);</span></div><div class='line' id='LC35'><span class="p">}</span></div><div class='line' id='LC36'><br/></div><div class='line' id='LC37'><span class="c1">//------------------------------------------------------------------------------</span></div><div class='line' id='LC38'><span class="n">vec</span> <span class="n">Blocking</span><span class="o">::</span><span class="n">readDataFromFile</span><span class="p">()</span></div><div class='line' id='LC39'><span class="p">{</span></div><div class='line' id='LC40'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vec</span> <span class="n">data</span><span class="p">;</span></div><div class='line' id='LC41'><br/></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Loading the first file</span></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">stringstream</span> <span class="n">file</span><span class="p">;</span></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">file</span> <span class="o">&lt;&lt;</span> <span class="n">dataPath</span> <span class="o">&lt;&lt;</span> <span class="s">&quot;/&quot;</span> <span class="o">&lt;&lt;</span> <span class="n">dataName</span> <span class="o">&lt;&lt;</span> <span class="mi">0</span> <span class="o">&lt;&lt;</span> <span class="s">&quot;.mat&quot;</span><span class="p">;</span></div><div class='line' id='LC45'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">data</span><span class="p">.</span><span class="n">load</span><span class="p">(</span><span class="n">file</span><span class="p">.</span><span class="n">str</span><span class="p">());</span></div><div class='line' id='LC46'><br/></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Loading the rest of the files</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vec</span> <span class="n">data_i</span><span class="p">;</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">nNodes</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">file</span><span class="p">.</span><span class="n">str</span><span class="p">(</span><span class="s">&quot;&quot;</span><span class="p">);</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">file</span> <span class="o">&lt;&lt;</span> <span class="n">dataPath</span> <span class="o">&lt;&lt;</span> <span class="s">&quot;/&quot;</span> <span class="o">&lt;&lt;</span> <span class="n">dataName</span> <span class="o">&lt;&lt;</span> <span class="n">i</span> <span class="o">&lt;&lt;</span> <span class="s">&quot;.mat&quot;</span><span class="p">;</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">data_i</span><span class="p">.</span><span class="n">load</span><span class="p">(</span> <span class="n">file</span><span class="p">.</span><span class="n">str</span><span class="p">()</span> <span class="p">);</span></div><div class='line' id='LC53'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">data</span> <span class="o">=</span> <span class="n">join_cols</span><span class="p">(</span><span class="n">data</span><span class="p">,</span> <span class="n">data_i</span><span class="p">);</span></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC55'><br/></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">data</span><span class="p">;</span></div><div class='line' id='LC57'><span class="p">}</span></div><div class='line' id='LC58'><span class="c1">//------------------------------------------------------------------------------</span></div><div class='line' id='LC59'><span class="n">vec</span> <span class="n">Blocking</span><span class="o">::</span><span class="n">computeBlock</span><span class="p">(</span><span class="k">const</span> <span class="n">vec</span> <span class="o">&amp;</span><span class="n">data</span><span class="p">,</span> <span class="kt">int</span> <span class="n">blockSize</span><span class="p">)</span></div><div class='line' id='LC60'><span class="p">{</span></div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">samples</span> <span class="o">=</span> <span class="n">data</span><span class="p">.</span><span class="n">size</span><span class="p">();</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">blocks</span> <span class="o">=</span> <span class="kt">int</span><span class="p">(</span><span class="n">samples</span> <span class="o">/</span> <span class="n">blockSize</span><span class="p">);</span></div><div class='line' id='LC63'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vec</span> <span class="n">averageEnergyBlock</span> <span class="o">=</span> <span class="n">zeros</span><span class="p">(</span><span class="n">blocks</span><span class="p">);</span></div><div class='line' id='LC64'><br/></div><div class='line' id='LC65'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Computing the mean value of every block</span></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">energySum</span><span class="p">;</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">blocks</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Finding the average over one block</span></div><div class='line' id='LC69'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">energySum</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC70'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">j</span> <span class="o">=</span> <span class="n">i</span> <span class="o">*</span> <span class="n">blockSize</span><span class="p">;</span> <span class="n">j</span> <span class="o">&lt;</span> <span class="n">i</span> <span class="o">*</span> <span class="n">blockSize</span> <span class="o">+</span> <span class="n">blockSize</span><span class="p">;</span> <span class="n">j</span><span class="o">++</span><span class="p">){</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">energySum</span> <span class="o">+=</span> <span class="n">data</span><span class="p">(</span><span class="n">j</span><span class="p">);</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">averageEnergyBlock</span><span class="p">(</span><span class="n">i</span><span class="p">)</span> <span class="o">=</span> <span class="n">energySum</span> <span class="o">/</span> <span class="n">blockSize</span><span class="p">;</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC75'><br/></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Calculating the mean and variance of all the blocks</span></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">E</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">E2</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC79'><br/></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">blocks</span><span class="p">;</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC81'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">E</span> <span class="o">+=</span> <span class="n">averageEnergyBlock</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">E2</span> <span class="o">+=</span> <span class="n">averageEnergyBlock</span><span class="p">(</span><span class="n">i</span><span class="p">)</span> <span class="o">*</span> <span class="n">averageEnergyBlock</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC83'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC84'><br/></div><div class='line' id='LC85'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Averaging</span></div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">E</span>  <span class="o">/=</span> <span class="n">blocks</span><span class="p">;</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">E2</span> <span class="o">/=</span> <span class="n">blocks</span><span class="p">;</span></div><div class='line' id='LC88'><br/></div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">sigma</span> <span class="o">=</span> <span class="n">E2</span> <span class="o">-</span> <span class="n">E</span><span class="o">*</span><span class="n">E</span><span class="p">;</span></div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="n">sigma</span><span class="o">/</span><span class="n">blocks</span><span class="p">);</span></div><div class='line' id='LC91'><br/></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vec</span> <span class="nf">result</span><span class="p">(</span><span class="mi">3</span><span class="p">);</span></div><div class='line' id='LC93'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">(</span><span class="mi">0</span><span class="p">)</span> <span class="o">=</span> <span class="n">blockSize</span><span class="p">;</span></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span> <span class="o">=</span> <span class="n">E</span><span class="p">;</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">(</span><span class="mi">2</span><span class="p">)</span> <span class="o">=</span> <span class="n">sigma</span><span class="p">;</span></div><div class='line' id='LC96'><br/></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">result</span><span class="p">;</span></div><div class='line' id='LC98'><span class="p">}</span></div><div class='line' id='LC99'><span class="c1">//------------------------------------------------------------------------------</span></div></pre></div>
          </td>
        </tr>
      </table>
  </div>

          </div>
        </div>

        <a href="#jump-to-line" rel="facebox" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
        <div id="jump-to-line" style="display:none">
          <h2>Jump to Line</h2>
          <form accept-charset="UTF-8" class="js-jump-to-line-form">
            <input class="textfield js-jump-to-line-field" type="text">
            <div class="full-button">
              <button type="submit" class="button">Go</button>
            </div>
          </form>
        </div>

      </div>
    </div>
</div>

<div id="js-frame-loading-template" class="frame frame-loading large-loading-area" style="display:none;">
  <img class="js-frame-loading-spinner" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-128.gif?1347543528" height="64" width="64">
</div>


        </div>
      </div>
      <div class="modal-backdrop"></div>
    </div>

      <div id="footer-push"></div><!-- hack for sticky footer -->
    </div><!-- end of wrapper - hack for sticky footer -->

      <!-- footer -->
      <div id="footer">
  <div class="container clearfix">

      <dl class="footer_nav">
        <dt>GitHub</dt>
        <dd><a href="/about">About us</a></dd>
        <dd><a href="/blog">Blog</a></dd>
        <dd><a href="/contact">Contact &amp; support</a></dd>
        <dd><a href="http://enterprise.github.com/">GitHub Enterprise</a></dd>
        <dd><a href="http://status.github.com/">Site status</a></dd>
      </dl>

      <dl class="footer_nav">
        <dt>Applications</dt>
        <dd><a href="http://mac.github.com/">GitHub for Mac</a></dd>
        <dd><a href="http://windows.github.com/">GitHub for Windows</a></dd>
        <dd><a href="http://eclipse.github.com/">GitHub for Eclipse</a></dd>
        <dd><a href="http://mobile.github.com/">GitHub mobile apps</a></dd>
      </dl>

      <dl class="footer_nav">
        <dt>Services</dt>
        <dd><a href="http://get.gaug.es/">Gauges: Web analytics</a></dd>
        <dd><a href="http://speakerdeck.com">Speaker Deck: Presentations</a></dd>
        <dd><a href="https://gist.github.com">Gist: Code snippets</a></dd>
        <dd><a href="http://jobs.github.com/">Job board</a></dd>
      </dl>

      <dl class="footer_nav">
        <dt>Documentation</dt>
        <dd><a href="http://help.github.com/">GitHub Help</a></dd>
        <dd><a href="http://developer.github.com/">Developer API</a></dd>
        <dd><a href="http://github.github.com/github-flavored-markdown/">GitHub Flavored Markdown</a></dd>
        <dd><a href="http://pages.github.com/">GitHub Pages</a></dd>
      </dl>

      <dl class="footer_nav">
        <dt>More</dt>
        <dd><a href="http://training.github.com/">Training</a></dd>
        <dd><a href="/edu">Students &amp; teachers</a></dd>
        <dd><a href="http://shop.github.com">The Shop</a></dd>
        <dd><a href="/plans">Plans &amp; pricing</a></dd>
        <dd><a href="http://octodex.github.com/">The Octodex</a></dd>
      </dl>

      <hr class="footer-divider">


    <p class="right">&copy; 2013 <span title="0.05233s from fe1.rs.github.com">GitHub</span>, Inc. All rights reserved.</p>
    <a class="left" href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>
    <ul id="legal">
        <li><a href="/site/terms">Terms of Service</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
    </ul>

  </div><!-- /.container -->

</div><!-- /.#footer -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/ComputationalPhysics/fys4411/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      Something went wrong with that request. Please try again.
      <a href="#" class="octicon octicon-remove-close ajax-error-dismiss"></a>
    </div>

    
    <span id='server_response_time' data-time='0.05271' data-host='fe1'></span>
    
  </body>
</html>

