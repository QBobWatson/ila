/*******************************************************************************
 * GTMathbook.js
 *******************************************************************************
 * The main front-end controller for GT Mathbook documents.
 *
 * Rewritten by Joseph Rabinoff.
 *
 * Original Authors: Michael DuBois, David Farmer, Rob Beezer
 *
 *******************************************************************************
 */

// Leading semicolon safeguards against errors in script concatenation
// Pass dependencies into this closure from the bottom of the file
;(function($, w) {
    'use strict'; // Use EMCAScript 5 strict mode within this closure

    // This class handles two things:
    //  (1) Resizing the nav dropdown and content relative to the viewport.
    //  (2) Toggling the nav dropdown.
    // In full-page mode, (2) shows a drop-down menu.  In mobile mode, (2)
    // replaces the page content with the menu; this allows for much easier
    // scrolling on mobile browsers.

    var Mathbook = function() {
        var $dropdown = $(".dropdown");
        var $tocContents = $(".toc-contents");
        var $tocBorder = $(".toc-border-container");
        var $toggle = $(".toggle-button");
        var sectionId = $("section").first().attr("id");
        var $sectionLink = $("a[data-scroll=" + sectionId + "]");
        var $topLink = $(".toc-contents h2.link").first();
        var $bottomButtons = $(".navbar-bottom-buttons");
        var $navbar = $("#gt-navbar");
        var $content = $("#content");
        var $main = $("main.main").first();
        var $w = $(w);

        var topNav = function() {
            return $bottomButtons.css("display") === "none";
        }

        var hideDropdown = function(e) {
            $tocBorder.hide();
            $toggle.removeClass("active");
            $content.show();
            if(topNav()) {
                $dropdown.append($tocBorder);
            } else {
                $main.prepend($tocBorder);
            }
        }

        var showDropdown = function(resizeOnly) {
            $toggle.addClass("active");
            if(topNav()) {
                $content.show();
                $dropdown.append($tocBorder);
            } else {
                $content.hide();
                $main.prepend($tocBorder);
            }
            $tocBorder.show();
            // Mobile browsers resize often based on which UI elements are
            // present.  We don't want to re-scroll every time.
            if(resizeOnly) {return}
            if(topNav()) {
                $tocContents.animate({
                    scrollTop: ($sectionLink.position().top
                                - $tocContents.height()/2
                                + $sectionLink.outerHeight(true)/2
                                - $topLink.position().top) + 'px'
                }, 'fast');
            } else {
                $("body").animate({
                    scrollTop: ($sectionLink.offset().top
                                + $sectionLink.outerHeight(true)
                                + $navbar.outerHeight(true)/2
                                - $w.height()/2) + 'px'
                }, 'fast');
            }
        }

        var toggleDropdown = function(e) {
            if($toggle.hasClass("active")) {
                hideDropdown();
            } else {
                showDropdown();
            }
        }

        var updateDropdown = function(resizeOnly) {
            if($toggle.hasClass("active")) {
                showDropdown(resizeOnly);
            } else {
                hideDropdown();
            }
        }

        var hideDropdownMaybe = function(e) {
            if(!e.target.matches(".dropdown") &&
               !e.target.matches(".toggle-button") &&
               $toggle.hasClass("active")) {
                hideDropdown();
            }
        }

        var tocHeight = function() {
            return ($w.height()
                    - $dropdown.offset().top + $(document).scrollTop()
                    - $tocBorder.outerHeight(true) + $tocBorder.innerHeight())
                * .85;
        }

        var resize = function(e) {
            // Stick or unstick the navbar based on it is at the top or bottom.
            // This avoids needing to reduntantly define the media sizes in JS.
            if(topNav()) {
                // Navbar at the top
                if(!$navbar.parent(".sticky-wrapper").length) {
                    $navbar.sticky();
                }
                $navbar.sticky("update");

                // Set the height of the ToC
                $tocContents.css({maxHeight: tocHeight()});
            } else {
                $navbar.unstick();
                $tocContents.css({maxHeight: 'none'});
                /*
                // Set the height of the ToC
                $tocContents.css({
                    maxHeight: $dropdown.offset().top
                        - $tocBorder.outerHeight(true)
                        + $tocBorder.innerHeight()
                });
                */
            }

            updateDropdown(true);

            // Make sure the content is large enough to vertically fill the
            // window.
            $content.css({
                minHeight: $w.height()
                    - ($("body").innerHeight() - $content.innerHeight())
            });
        }

        var onScroll = function(e) {
            if(topNav()) {
                // Set the height of the ToC (which may have changed)
                $tocContents.css({maxHeight: tocHeight()});
            }
        }

        $toggle.on("click", toggleDropdown);
        $w.on("click", hideDropdownMaybe);
        $w.on("resize", resize);
        $w.scroll(onScroll);
        resize();

        // Hack
        if(!$("#toc a.active").length) {
            $("#toc h2.active a").addClass("active");
        }

        $(".mathbook-content section.hidden-subsection > header > h1").on(
            'click', function() {
                var parent = $(this).parent().parent();
                var child = parent.children(".hidden-subsection-content");
                if(parent.hasClass("active")) {
                    parent.removeClass("active");
                } else {
                    parent.addClass("active");
                }
                child.slideToggle(500);
            });

        var maximize = function() {
            var parent = $(this).parent();
            var iframe = parent.children("iframe");
            if(parent.hasClass("maximized")) {
                parent.removeClass("maximized");
                parent.css("height", parent.data("height"));

                var iOS = !!navigator.platform &&
                    /iPad|iPhone|iPod/.test(navigator.platform);
                if(iOS) {
                    // Hack for iOS
                    iframe.remove();
                    parent.prepend(iframe);
                }
            } else {
                parent.data("height", parent.css("height"));
                parent.css("height", "auto");
                parent.addClass("maximized");
            }
        };

        $(".mathbook-content .mathbox .maximizer")
            .on("click", maximize);
        $(".mathbook-content .mathbox .minimizer")
            .on("click", maximize);

        this.knowlLoaded = function(uid) {
            $("#kuid-" + uid + " .mathbox .maximizer")
                .on("click", maximize);
            $("#kuid-" + uid + " .mathbox .minimizer")
                .on("click", maximize);
            $("#kuid-" + uid + " .concept-button")
                .click(function() {
                    $(this).parent().toggleClass("show");
                });
        };

        // Concept library
        $(".concept-button").click(function() {
            $(this).parent().toggleClass("show");
        });
    };

    // If script is run after page is loaded, initialize immediately
    if(document.readyState === "complete") {
        w.mathbook = new Mathbook();
    } else {
        // wait and init when the DOM is fully loaded
        $(window).load( function() {
            w.mathbook = new Mathbook();
        });
    }

    return Mathbook;

})(jQuery, window);
