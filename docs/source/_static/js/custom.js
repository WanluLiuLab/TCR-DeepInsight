document.querySelectorAll('.sig-paren').forEach(el => {
    if (el.textContent.trim() === ')') {
        if ((hasEmInHierarchy = $(el).parent().find("em").length > 0) & ($(el).parent().text().includes(', '))) {
            el.style.setProperty('--line-break-content', 'block');
        }
    }
});

$(".theme-icon-when-auto-dark").remove();