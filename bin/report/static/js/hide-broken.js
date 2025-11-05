function hideBrokenLinks() {
  // Find all elements with class 'hide-broken' and hide any child elements with
  // an href that returns an error code
  const elements = document.getElementsByClassName('hide-broken');
  Array.from(elements).forEach(element => {
    // Get all child elements with an href attribute
    const children = element.querySelectorAll('[href]');
    const hrefElements = [element].concat(Array.from(children));
    hrefElements.forEach(el => {
      // Check if the href returns an error code
      if (!el.href) {
        return;
      }
      hideIfLinkBroken(element, el.href);
    });
  });
}

function hideIfLinkBroken(element, url) {
  // Unfortunately, this can't be determined with an offline file:// HTML page
  return false;
}

document.addEventListener('DOMContentLoaded', hideBrokenLinks);
