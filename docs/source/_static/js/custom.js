document.querySelectorAll('.sig-paren').forEach(el => {
    if (el.textContent.trim() === ')') {
        if (Array.from(
                $(el).parent().find("em")
        ).filter(e => !e.className.includes("property")).length > 1) {
            el.style.setProperty('--line-break-content', 'block');
        }
    }
});
document.querySelectorAll('.sig-param').forEach(el => {
    if (Array.from($(el).parent().find("em")).filter(e => !e.className.includes("property")).length < 2) {
        $(el).addClass("remove-before")
    } else {
        var t = $($(el).find(".n")[0]).text();
        $($(el).find(".n")[0]).text("    " + t)
    }
})


$(".theme-icon-when-auto-dark").remove();