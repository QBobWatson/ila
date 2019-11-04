/*
 * Knowl - Feature Demo for Knowls
 * Copyright (C) 2011  Harald Schilly
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * 4/11/2012 Modified by David Guichard to allow inline knowl code.
 * Sample use:
 *      This is an <a knowl="" class="internal"
 *      value="Hello World!">inline knowl.</a>
 */

/*  8/14/14  Modified by David Farmer to allow knowl content to be
 *  taken from the element with a given id.
 *
 * The syntax is <a knowl="" class="id-ref" refid="proofSS">Proof</a>
 */

/* javascript code for the knowl features
 * global counter, used to uniquely identify each knowl-output element
 * that's necessary because the same knowl could be referenced several times
 * on the same page */
var knowl_id_counter = 0;

var knowl_focus_stack_uid = [];
var knowl_focus_stack = [];

function knowl_click_handler($el) {
  // the knowl attribute holds the id of the knowl
  var knowl_id = $el.attr("knowl");
  // the uid is necessary if we want to reference the same content several times
  var uid = $el.attr("knowl-uid");
  var output_id = '#knowl-output-' + uid;
  var $output_id = $(output_id);
  // create the element for the content, insert it after the one where the
  // knowl element is included (e.g. inside a <h1> tag) (sibling in DOM)
  var idtag = "id='"+output_id.substring(1) + "'";
  var kid   = "id='kuid-"+ uid + "'";
  // if we already have the content, toggle visibility

  // Note that for tracking knowls, this setup is not optimal
  // because it applies to open knowls and also knowls which
  // were opened and then closed.
  if ($output_id.length > 0) {
     thisknowlid = "kuid-"+uid
     $("#kuid-"+uid).slideToggle("fast");
     if($el.attr("replace")) {
       $($el.attr("replace")).slideToggle("fast");
     }

     this_knowl_focus_stack_uidindex = knowl_focus_stack_uid.indexOf(uid);

     if($el.hasClass("active")) {
       if(this_knowl_focus_stack_uidindex != -1) {
         knowl_focus_stack_uid.splice(this_knowl_focus_stack_uidindex, 1);
         knowl_focus_stack.splice(this_knowl_focus_stack_uidindex, 1);
       }
     }
     else {
         knowl_focus_stack_uid.push(uid);
         knowl_focus_stack.push($el);
         document.getElementById(thisknowlid).focus();
     }

     $el.toggleClass("active");

  // otherwise download it or get it from the cache
  } else {
    // where_it_goes is the location the knowl will appear *after*
    // knowl is the variable that will hold the content of the output knowl
    var where_it_goes = $el;
    var knowl = "<div class='knowl-output' "+kid+"><div class='knowl'><div class='knowl-content' " +idtag+ ">loading '"+knowl_id+"'</div><div class='knowl-footer'>"+knowl_id+"</div></div></div>";

    // addafter="#id" means to put the knowl after the element with that id
    if($el.attr("addafter")) {
        where_it_goes = $($el.attr("addafter"));
    } else if($el.attr("replace")) {
        where_it_goes = $($el.attr("replace"));
    } else if($el.hasClass("kohere")) {
        where_it_goes = $el;
    } else {
       // otherwise, typically put it after the nearest enclosing block element

      // check, if the knowl is inside a td or th in a table
      if($el.parent().is("td") || $el.parent().is("th") ) {
        // assume we are in a td or th tag, go 2 levels up
        where_it_goes = $el.parent().parent();
        var cols = $el.parent().parent().children().length;
        knowl = "<tr><td colspan='"+cols+"'>"+knowl+"</td></tr>";
      } else if ($el.parent().is("li")) {
        where_it_goes = $el.parent();
      }
      // not sure it is is worth making the following more elegant
      else if ($el.parent().parent().is("li")) {
        where_it_goes = $el.parent().parent();
        // the '.is("p")' is for the first paragraph of a theorem or proof
      } else if ($el.parent().css('display') == "block" || $el.parent().is("p") || $el.parent().hasClass("hidden-knowl-wrapper") || $el.parent().hasClass("kohere")) {
        where_it_goes = $el.parent();
      } else if ($el.parent().parent().css('display') == "block" || $el.parent().parent().is("p") || $el.parent().parent().hasClass("hidden-knowl-wrapper") || $el.parent().parent().hasClass("kohere")) {
        where_it_goes = $el.parent().parent();
      } else {
        //  is this a reasonable last case?
        //  if we omit the else, then if goes after $el
        where_it_goes = $el.parent().parent().parent();
      }

    }

    // now that we know where the knowl goes, insert the knowl content
    if($el.attr("replace")) {
        where_it_goes.before(knowl);
    }
    else {
        where_it_goes.after(knowl);
    }

    // "select" where the output is and get a hold of it
    var $output = $(output_id);
    var $knowl = $("#kuid-"+uid);
    $output.addClass("loading");
    $knowl.hide();

    // DRG: inline code
    if ($el.hasClass('internal')) {
      $output.html($el.attr("value"));
//    } else if ($el.attr("class") == 'id-ref') {
    } else if ($el.hasClass('id-ref')) {
     //get content from element with the given id
      $output.html($("#".concat($el.attr("refid"))).html());
    } else {
    // Get code from server.
    $output.load(knowl_id,
     function(response, status, xhr) {
       $knowl.removeClass("loading");
       if (status == "error") {
         $el.removeClass("active");
         $output.html("<div class='knowl-output error'>ERROR: " + xhr.status + " " + xhr.statusText + '</div>');
         $output.show();
       } else if (status == "timeout") {
         $el.removeClass("active");
         $output.html("<div class='knowl-output error'>ERROR: timeout. " + xhr.status + " " + xhr.statusText + '</div>');
         $output.show();
       } else {
         // Success
         window.mathbook.knowlLoaded(uid);
    }
     });
    };

   // we have the knowl content, and put it hidden in the right place,
   // so now we show it

   $knowl.hide();

   $el.addClass("active");
   $knowl.slideDown("slow");
  }
} //~~ end click handler for *[knowl] elements

/** register a click handler for each element with the knowl attribute
 * @see jquery's doc about 'live'! the handler function does the
 *  download/show/hide magic. also add a unique ID,
 *  necessary when the same reference is used several times. */
$(function() {
    $("body").on("click", "*[knowl]", function(evt) {
      evt.preventDefault();
      var $knowl = $(this);
      if(!$knowl.attr("knowl-uid")) {
        $knowl.attr("knowl-uid", knowl_id_counter);
        knowl_id_counter++;
      }
      knowl_click_handler($knowl, evt);
  });
});


$(window).load(function() {
   $("a[knowl]").attr("href", "");
});

//window.onload = function() {
/*
window.addEventListener("load",function(event) {
    document.onkeyup = function(event)
    {
        var e = (!event) ? window.event : event;
        switch(e.keyCode)
        {
            case 27: //esc
                if(knowl_focus_stack.length > 0 ) {
                   most_recently_opened = knowl_focus_stack.pop();
                   knowl_focus_stack_uid.pop();
                   most_recently_opened.focus();
                } else {
                   console.log("no open knowls being tracked");
                   break;
                }
        };
    };
},
false);

*/
