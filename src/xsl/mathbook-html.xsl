<?xml version='1.0'?>  <!-- -*- nxml-child-indent: 4 -*- -->

<!DOCTYPE xsl:stylesheet [
    <!ENTITY % entities SYSTEM "../../mathbook/xsl/entities.ent">
    %entities;
]>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
    xmlns:exsl="http://exslt.org/common"
    extension-element-prefixes="exsl"
>

<xsl:import href="../../mathbook/xsl/mathbook-html.xsl" />
<xsl:import href="git-hash.xsl" />
<xsl:import href="theme.xsl" />

<!-- JDR: for caching -->
<xsl:param name="debug.datedfiles">no</xsl:param>

<xsl:param name="toc.level" select="2" />
<!-- <xsl:param name="html.knowl.example" select="'no'" /> -->

<xsl:param name="pdf.online" select="'ila.pdf'" />

<!-- JDR: one js file to rule them all -->
<xsl:template name="combined-js">
    <xsl:element name="script">
        <xsl:attribute name="src">
            <xsl:call-template name="versioned-file">
                <xsl:with-param name="file" select="'js/ila.js'"/>
            </xsl:call-template>
        </xsl:attribute>
    </xsl:element>
</xsl:template>

<xsl:template name="mathbook-js"/>
<xsl:template name="jquery-sagecell"/>

<!-- JDR: one css file to rule them all -->
<xsl:template name="css">
    <xsl:element name="link">
        <xsl:attribute name="href">
            <xsl:call-template name="versioned-file">
                <xsl:with-param name="file" select="'css/ila.css'"/>
            </xsl:call-template>
        </xsl:attribute>
        <xsl:attribute name="rel">stylesheet</xsl:attribute>
        <xsl:attribute name="type">text/css</xsl:attribute>
    </xsl:element>
    <xsl:element name="link">
        <xsl:attribute name="rel">manifest</xsl:attribute>
        <xsl:attribute name="href">manifest.json</xsl:attribute>
    </xsl:element>
    <!-- JDR: preprocessed inline pretex stylesheet is inserted here -->
    <style id="pretex-style"></style>
    <style id="pretex-fonts"></style>
    <link rel="shortcut icon" href="images/logo.png"/>
</xsl:template>

<!-- JDR: use theme.sponsor in brand logo -->
<xsl:template name="brand-logo">
    <xsl:choose>
        <xsl:when test="/mathbook/docinfo/brandlogo">
            <a id="logo-link" href="{$theme.sponsor}" target="_blank" rel="noopener">
                <img src="{/mathbook/docinfo/brandlogo/@source}" alt="Logo image for document"/>
            </a>
        </xsl:when>
        <xsl:otherwise>
            <a id="logo-link" href="" />
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<!-- Primary Navigation -->
<!-- ToC, Prev/Up/Next buttons  -->
<xsl:template match="*" mode="primary-navigation">
  <nav id="gt-navbar" class="navbar container" style="">
    <div class="dropdown">
      <div class="toc-border-container" id="toc">
        <div class="toc-contents">
          <xsl:apply-templates select="." mode="toc-items" />
        </div>
      </div>
    </div>
    <div class="navbar-top-buttons toolbar">
      <div class="toolbar-buttons-left">
        <!-- Toggle button -->
        <a class="toggle-button button toolbar-item"
           title="Toggle Menu" alt="Toggle Menu"
           href="javascript:void(0)"></a>
        <!-- A page either has an/the index as    -->
        <!-- a child, and gets the "jump to" bar, -->
        <!-- or it deserves an index button       -->
        <xsl:choose>
          <xsl:when test="index-list">
            <div class="toolbar-item">
              <xsl:apply-templates select="." mode="index-jump-nav" />
            </div>
          </xsl:when>
          <xsl:otherwise>
            <xsl:apply-templates select="." mode="index-button" />
          </xsl:otherwise>
        </xsl:choose>
      </div>
      <div class="toolbar-buttons-right">
        <!-- Each button gets an id for keypress recognition/action -->
        <xsl:apply-templates select="." mode="previous-button">
          <xsl:with-param name="id-label" select="'previousbutton-top'" />
        </xsl:apply-templates>
        <xsl:if test="$nav-upbutton='yes'">
          <xsl:apply-templates select="." mode="up-button">
            <xsl:with-param name="id-label" select="'upbutton-top'" />
          </xsl:apply-templates>
        </xsl:if>
        <xsl:apply-templates select="." mode="next-button">
          <xsl:with-param name="id-label" select="'nextbutton-top'" />
        </xsl:apply-templates>
      </div>
    </div>
    <div class="navbar-bottom-buttons toolbar">
      <a class="toggle-button button toolbar-item"
         title="Toggle Menu" alt="Toggle Menu"
         href="javascript:void(0)"></a>
      <xsl:apply-templates select="." mode="previous-button">
        <xsl:with-param name="id-label" select="'previousbutton-bot'" />
      </xsl:apply-templates>
      <xsl:if test="$nav-upbutton='yes'">
        <xsl:apply-templates select="." mode="up-button">
          <xsl:with-param name="id-label" select="'upbutton-bot'" />
        </xsl:apply-templates>
      </xsl:if>
      <xsl:apply-templates select="." mode="next-button">
        <xsl:with-param name="id-label" select="'nextbutton-bot'" />
      </xsl:apply-templates>
    </div>
  </nav>
</xsl:template>

<!-- JDR: no sidebars anymore -->
<xsl:template match="*" mode="sidebars">
</xsl:template>

<!-- JDR: "answers" in examples not hidden -->
<xsl:template match="example/hint|example/answer|example/solution" mode="is-hidden">
    <xsl:text>false</xsl:text>
</xsl:template>
<xsl:template match="specialcase/hint|specialcase/answer|specialcase/solution" mode="is-hidden">
    <xsl:text>false</xsl:text>
</xsl:template>

<!-- JDR: un-hide a proof -->
<xsl:template match="proof" mode="is-hidden">
    <xsl:choose>
        <xsl:when test="@visible='true'">
            <xsl:value-of select="false()" />
        </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="$html.knowl.proof = 'yes'" />
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<!-- JDR: "special case" is just a non-hidden example -->
<xsl:template match="specialcase" mode="is-hidden">
    <xsl:text>false</xsl:text>
</xsl:template>

<!-- JDR: remarks are hidden -->
<xsl:template match="remark" mode="is-hidden">
    <xsl:text>true</xsl:text>
</xsl:template>

<!-- JDR: paragraphs are optionally hidden -->
<xsl:template match="paragraphs" mode="is-hidden">
    <xsl:choose>
        <xsl:when test="@visible='false'">
            <xsl:text>true</xsl:text>
        </xsl:when>
        <xsl:otherwise>
            <xsl:text>false</xsl:text>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template match="bluebox" mode="body-css-class">
    <xsl:text>bluebox</xsl:text>
</xsl:template>

<!-- JDR: this produces the links in the notation list and concept library -->

<xsl:template match="*" mode="enclosure-xref">
    <xsl:variable name="structural">
        <xsl:apply-templates select="." mode="is-structural" />
    </xsl:variable>
    <xsl:variable name="block">
        <xsl:apply-templates select="." mode="is-block" />
    </xsl:variable>
    <xsl:choose>
        <!-- found a structural or block parent -->
        <!-- we fashion a cross-reference link  -->
        <xsl:when test="$structural='true' or $block='true'">
            <xsl:apply-templates select="." mode="xref-link">
                <xsl:with-param name="content">
                    <xsl:apply-templates select="." mode="type-name" />
                </xsl:with-param>
            </xsl:apply-templates>
        </xsl:when>
        <!-- nothing interesting here, so step up a level -->
        <!-- Eventually we find the top-level structure   -->
        <!-- eg article, book, etc                        -->
        <xsl:otherwise>
            <xsl:apply-templates select="parent::*" mode="enclosure-xref" />
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<!-- JDR: Concept library support -->

<xsl:template match="concept" />

<xsl:template match="*" mode="concept">
    <xsl:if test="concept">
        <span class="concepts">
            <span class="concept-button">&#xf2e8;</span>
            <span class="concept-list">
                <xsl:for-each select="concept">
                    <xsl:variable name="tag">
                        <xsl:value-of select="@tag" />
                    </xsl:variable>
                    <xsl:if test="not(//concept-library/concept[@tag=$tag])">
                        <xsl:message>MBX:WARNING: Concept library tag #<xsl:value-of select="@tag" /> is used but not defined in the concept library.</xsl:message>
                    </xsl:if>
                    <xsl:element name="a">
                        <xsl:attribute name="href">
                            <xsl:text>concept-library.html#concept-</xsl:text>
                            <xsl:value-of select="@tag" />
                        </xsl:attribute>
                        <xsl:text>#</xsl:text><xsl:value-of select="@tag" />
                    </xsl:element>
                    <xsl:if test="position() != last()">
                        <br/>
                    </xsl:if>
                </xsl:for-each>
            </span>
        </span>
    </xsl:if>
</xsl:template>

<xsl:template match="concept-library">
    <table class="concept-library">
        <xsl:apply-templates select="concept" mode="library" />
    </table>
</xsl:template>

<xsl:template match="concept" mode="library">
    <xsl:variable name="tag">
        <xsl:value-of select="@tag"/>
    </xsl:variable>
    <tr>
        <xsl:element name="td">
            <xsl:attribute name="class">
                <xsl:text>concept-library-tag</xsl:text>
            </xsl:attribute>
            <xsl:attribute name="id">
                <xsl:text>concept-</xsl:text><xsl:value-of select="@tag"/>
            </xsl:attribute>
            <xsl:text>#</xsl:text><xsl:value-of select="@tag"/>
        </xsl:element>
        <td>
            <xsl:apply-templates select="desc" />
        </td>
        <td>
            <xsl:apply-templates select="//concept[@tag=$tag][not(parent::concept-library)]" mode="enclosure-xref" />
            <xsl:if test="see">
                <br/>
                <i class="see">also see</i>
                <xsl:for-each select="see">
                    <xsl:variable name="see-tag">
                        <xsl:value-of select="@tag" />
                    </xsl:variable>
                    <xsl:if test="not(//concept-library/concept[@tag=$see-tag])">
                        <xsl:message>MBX:WARNING: Concept library tag #<xsl:value-of select="@tag" /> is used but not defined in the concept library.</xsl:message>
                    </xsl:if>
                    <xsl:element name="a">
                        <xsl:attribute name="class">
                            <xsl:text>concept-library-tag</xsl:text>
                        </xsl:attribute>
                        <xsl:attribute name="href">
                            <xsl:text>#concept-</xsl:text>
                            <xsl:value-of select="@tag" />
                        </xsl:attribute>
                        <xsl:text>#</xsl:text><xsl:value-of select="@tag" />
                    </xsl:element>
                    <xsl:if test="position() != last()">
                        <xsl:text>, </xsl:text>
                    </xsl:if>
                </xsl:for-each>
            </xsl:if>
        </td>
    </tr>
</xsl:template>


<!-- JDR: mathbox support -->

<xsl:template match="mathbox" mode="panel-html-box">
    <xsl:apply-templates select="." />
</xsl:template>

<xsl:template match="mathbox">
    <xsl:variable name="height">
      <xsl:choose>
        <xsl:when test="@height">
          <xsl:value-of select="@height"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>300px</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="src">
      <xsl:choose>
        <xsl:when test="contains(@source, '?')">
          <xsl:value-of select="@source"/>
          <xsl:text>&amp;vers=</xsl:text>
          <xsl:call-template name="git-hash-short" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="@source"/>
          <xsl:text>?vers=</xsl:text>
          <xsl:call-template name="git-hash-short" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:element name="div">
        <xsl:attribute name="id">
            <xsl:apply-templates select="." mode="internal-id" />
        </xsl:attribute>
        <xsl:attribute name="class">
            <xsl:text>mathbox</xsl:text>
        </xsl:attribute>
        <xsl:attribute name="style">
            <xsl:text>height:</xsl:text>
            <xsl:value-of select="$height" />
            <xsl:text>;</xsl:text>
        </xsl:attribute>

        <xsl:element name="iframe">
          <xsl:attribute name="src">
            <xsl:value-of select="$src"/>
          </xsl:attribute>
          <xsl:attribute name="title">Interactive Figure</xsl:attribute>
        </xsl:element>

        <xsl:element name="div">
          <xsl:attribute name="class">
            <xsl:text>maximizer</xsl:text>
          </xsl:attribute>
          <xsl:text>&#xf25e;</xsl:text>
        </xsl:element>
        <xsl:element name="div">
          <xsl:attribute name="class">
            <xsl:text>minimizer</xsl:text>
          </xsl:attribute>
          <xsl:text>&#xf267;</xsl:text>
        </xsl:element>
    </xsl:element>

    <!--
    <xsl:element name="div">
      <xsl:attribute name="class">
        <xsl:text>mathbox-link</xsl:text>
      </xsl:attribute>
      <xsl:element name="a">
        <xsl:attribute name="href">
          <xsl:value-of select="$src"/>
        </xsl:attribute>
        <xsl:attribute name="target">_blank</xsl:attribute>
        <xsl:attribute name="rel">noopener</xsl:attribute>
        <xsl:text>Click to view in a new window</xsl:text>
      </xsl:element>
    </xsl:element>
    -->
</xsl:template>

<xsl:template match="*" mode="get-hide-type">
    <xsl:value-of select="@hide-type"/>
</xsl:template>

<!-- Blue boxes hide type by default -->
<xsl:template match="bluebox" mode="get-hide-type">
    <xsl:choose>
        <xsl:when test="@hide-type">
            <xsl:value-of select="@hide-type"/>
        </xsl:when>
        <xsl:otherwise>
            <xsl:text>true</xsl:text>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>


<!-- JDR: simpler numbering of some elements -->
<xsl:template match="*" mode="heading-simple-nonumber">
    <xsl:param name="important"/>
    <xsl:param name="in-knowl"/>
    <xsl:variable name="hide-type">
        <xsl:apply-templates select="." mode="get-hide-type"/>
    </xsl:variable>
    <xsl:if test="title or $hide-type != 'true'">
        <xsl:element name="h5">
            <xsl:attribute name="class">
                <xsl:text>heading</xsl:text>
                <xsl:if test="$important">
                   <xsl:text> important</xsl:text>
                </xsl:if>
            </xsl:attribute>
            <xsl:if test="$important">
                <xsl:element name="img">
                    <xsl:attribute name="src">
                        <xsl:text>images/important.svg</xsl:text>
                    </xsl:attribute>
                    <xsl:attribute name="class">
                        <xsl:text>important</xsl:text>
                    </xsl:attribute>
                    <xsl:attribute name="alt">Exclamation Point</xsl:attribute>
                </xsl:element>
            </xsl:if>
            <xsl:choose>
                <xsl:when test="$hide-type = 'true'">
                    <span class="title">
                        <xsl:apply-templates select="." mode="title-full" />
                    </span>
                </xsl:when>
                <xsl:otherwise>
                    <span class="type">
                        <xsl:apply-templates select="." mode="type-name" />
                    </span>
                    <xsl:if test="title">
                        <span class="title">
                            <xsl:text>(</xsl:text>
                            <xsl:apply-templates select="." mode="title-full" />
                            <xsl:text>)</xsl:text>
                        </span>
                    </xsl:if>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:element>
    </xsl:if>
</xsl:template>

<!-- JDR: customize heading-full -->
<xsl:template match="*" mode="heading-full">
    <xsl:param name="important"/>
    <xsl:element name="h5">
        <xsl:attribute name="class">
            <xsl:text>heading</xsl:text>
            <xsl:if test="$important">
               <xsl:text> important</xsl:text>
            </xsl:if>
        </xsl:attribute>
        <xsl:if test="$important">
            <xsl:element name="img">
                <xsl:attribute name="src">
                    <xsl:text>images/important.svg</xsl:text>
                </xsl:attribute>
                <xsl:attribute name="class">
                    <xsl:text>important</xsl:text>
                </xsl:attribute>
                <xsl:attribute name="alt">Exclamation Point</xsl:attribute>
            </xsl:element>
        </xsl:if>
        <xsl:variable name="the-number">
            <xsl:apply-templates select="." mode="number" />
        </xsl:variable>
        <span class="type">
            <xsl:apply-templates select="." mode="type-name" />
        </span>
        <xsl:if test="not($the-number='')">
            <span class="codenumber">
                <xsl:value-of select="$the-number" />
            </span>
        </xsl:if>
        <xsl:if test="title">
            <span class="title">
                <xsl:text>(</xsl:text>
                <xsl:apply-templates select="." mode="title-full" />
                <xsl:text>)</xsl:text>
            </span>
        </xsl:if>
    </xsl:element>
</xsl:template>

<xsl:template match="&DEFINITION-LIKE;|&REMARK-LIKE;" mode="heading-birth">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="&DEFINITION-LIKE;|&REMARK-LIKE;" mode="heading-xref-knowl">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="essential" mode="heading-birth">
    <xsl:apply-templates select="." mode="heading-simple-nonumber">
        <xsl:with-param name="important" select="true()"/>
    </xsl:apply-templates>
</xsl:template>

<xsl:template match="essential" mode="heading-xref-knowl">
    <xsl:apply-templates select="." mode="heading-simple-nonumber">
        <xsl:with-param name="important" select="true()"/>
    </xsl:apply-templates>
</xsl:template>

<xsl:template match="&EXAMPLE-LIKE;|&PROJECT-LIKE;|list" mode="heading-birth">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="&EXAMPLE-LIKE;|&PROJECT-LIKE;|list" mode="heading-xref-knowl">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="&THEOREM-LIKE;|&AXIOM-LIKE;" mode="heading-birth">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="&THEOREM-LIKE;|&AXIOM-LIKE;" mode="heading-xref-knowl">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>

<xsl:template match="exercise" mode="heading-birth">
    <xsl:apply-templates select="." mode="heading-simple-nonumber" />
</xsl:template>


<!-- JDR: do not preload example knowls.                       -->
<!-- This prevents initializing all the mathboxes on the page. -->
<xsl:template match="example" mode="born-hidden">
    <xsl:param name="b-original" select="true()" />
    <xsl:variable name="birth-elt">
        <xsl:apply-templates select="." mode="birth-element" />
    </xsl:variable>
    <!-- First: the link that is visible on the page -->
    <xsl:element name="{$birth-elt}">
        <xsl:attribute name="class">
            <xsl:text>hidden-knowl-wrapper</xsl:text>
        </xsl:attribute>
        <xsl:apply-templates select="." mode="concept" />
        <xsl:element name="a">
            <xsl:attribute name="knowl">
                <xsl:apply-templates select="." mode="hidden-knowl-url" />
            </xsl:attribute>
            <xsl:attribute name="knowl-id">
                <xsl:text>hidden-</xsl:text>
                <xsl:apply-templates select="." mode="internal-id" />
            </xsl:attribute>
            <!-- add HTML title and alt attributes to the link -->
            <xsl:attribute name="alt">
                <xsl:apply-templates select="." mode="tooltip-text" />
            </xsl:attribute>
            <xsl:attribute name="title">
                <xsl:apply-templates select="." mode="tooltip-text" />
            </xsl:attribute>
            <xsl:apply-templates select="." mode="hidden-knowl-text" />
        </xsl:element>
    </xsl:element>
</xsl:template>

<xsl:template name="knowl">
</xsl:template>

<!-- JDR: we're using CharterBT -->
<xsl:template name="fonts">
</xsl:template>

<!-- JDR: we're precompiling latex -->
<xsl:template name="mathjax"/>
<xsl:template name="latex-macros">
  <exsl:document href="./preamble.tex" method="text">
    <xsl:value-of select="$latex-packages-mathjax" />
    <xsl:value-of select="$latex-macros" />
  </exsl:document>
</xsl:template>

<!-- JDR: pretex can compile (almost) arbitrary latex code! -->
<!-- This is almost the same as <me>, except it doesn't wrap the result in an
     equation* environtment.  Also, in "bare" mode, commands can be executed in
     a global context. -->
<xsl:template match="latex-code">
  <xsl:choose>
    <xsl:when test="@mode='inline'">
      <script type="text/x-latex-code-inline">
        <xsl:value-of select="text()" />
      </script>
    </xsl:when>
    <xsl:when test="@mode='bare'">
      <script type="text/x-latex-code-bare">
        <xsl:value-of select="text()" />
      </script>
    </xsl:when>
    <xsl:otherwise>
      <script type="text/x-latex-code">
        <xsl:value-of select="text()" />
      </script>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- JDR: support for hidden (not knowl-ized!) subsections. -->
<xsl:template match="subsection">
    <!-- location info for debugging efforts -->
    <xsl:apply-templates select="." mode="debug-location" />
    <!-- Heading, div for this structural subdivision -->
    <xsl:variable name="ident">
      <xsl:apply-templates select="." mode="internal-id" />
    </xsl:variable>
    <xsl:choose>
      <xsl:when test="@hidden='true'">
        <section class="{local-name(.)} hidden-subsection" id="{$ident}">
          <xsl:apply-templates select="." mode="section-header" />
          <div class="hidden-subsection-content">
            <xsl:apply-templates />
          </div>
        </section>
      </xsl:when>
      <xsl:otherwise>
        <section class="{local-name(.)}" id="{$ident}">
          <xsl:apply-templates select="." mode="section-header" />
          <xsl:apply-templates />
        </section>
      </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<!-- JDR: feedback link -->
<xsl:template name="feedback-link">
  <div class="feedback-link">
    <xsl:element name="a">
      <xsl:attribute name="href">
        <xsl:text>https://github.com/QBobWatson/ila/issues/new?title=Feedback%20on%20</xsl:text>
        <xsl:apply-templates select="." mode="internal-id" />
        <xsl:text>&amp;body=In%20reference%20to%20version%20</xsl:text>
        <xsl:call-template name="git-hash" />
      </xsl:attribute>
      <xsl:attribute name="class">
        <xsl:text>feedback-link</xsl:text>
      </xsl:attribute>
      <xsl:attribute name="target">_blank</xsl:attribute>
      <xsl:attribute name="rel">noopener</xsl:attribute>
      <xsl:text>Comments, corrections or suggestions?</xsl:text>
      <span class="feedback-link">(Free GitHub account required)</span>
    </xsl:element>
  </div>
</xsl:template>

<!-- JDR: pdf version -->
<xsl:template name="pdf-version">
  <div class="pdf-version">
    <xsl:element name="a">
      <xsl:attribute name="href">
        <xsl:value-of select="$pdf.online" />
      </xsl:attribute>
      <xsl:attribute name="target">_blank</xsl:attribute>
      <xsl:attribute name="rel">noopener</xsl:attribute>
      <xsl:text>PDF version</xsl:text>
    </xsl:element>
  </div>
</xsl:template>

<!-- JDR: number overrides -->
<xsl:template match="part[@number]" mode="raw-serial-number">
    <xsl:number format="I" value="@number" />
</xsl:template>
<!-- TODO: condition on part/chapter style to use  level='any'; from="book/part"  to cross part boundaries -->
<xsl:template match="chapter[@number]" mode="raw-serial-number">
    <xsl:number value="@number" format="1" />
</xsl:template>
<xsl:template match="appendix[@number]" mode="raw-serial-number">
    <xsl:number value="@number" format="A" />
</xsl:template>
<xsl:template match="section[@number]" mode="raw-serial-number">
    <xsl:number value="@number" format="1" />
</xsl:template>
<xsl:template match="subsection[@number]" mode="raw-serial-number">
    <xsl:number value="@number" format="1" />
</xsl:template>
<xsl:template match="subsubsection[@number]" mode="raw-serial-number">
    <xsl:number value="@number" format="1" />
</xsl:template>


<!-- JDR: cover animation -->
<xsl:template name="cover">
  <div id="cover">
    <xsl:element name="script">
        <xsl:attribute name="src">
            <xsl:call-template name="versioned-file">
                <xsl:with-param name="file" select="'js/cover.js'"/>
            </xsl:call-template>
        </xsl:attribute>
    </xsl:element>
  </div>
</xsl:template>

</xsl:stylesheet>
