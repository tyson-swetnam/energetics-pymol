/* theme.js — light / dark / system theme switch, shared by index.html and learn.html.
   The actual colours come from light-dark() in style.css; this just sets color-scheme
   (via data-theme on <html>), persists the choice, and updates the toggle + any viewer. */
(function () {
  var html = document.documentElement;

  function resolvedDark(t) {
    return t === "dark" ||
      (t === "system" && window.matchMedia &&
       window.matchMedia("(prefers-color-scheme: dark)").matches);
  }

  function apply(t) {
    html.setAttribute("data-theme", t);
    html.style.colorScheme = (t === "system" ? "light dark" : t);
    try { localStorage.setItem("theme", t); } catch (e) {}
    document.querySelectorAll(".theme-toggle button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.dataset.themeSet === t));
    });
    if (typeof window.onThemeChange === "function") window.onThemeChange(resolvedDark(t));
  }

  var current = (function () {
    try { return localStorage.getItem("theme"); } catch (e) { return null; }
  })() || "system";

  document.querySelectorAll(".theme-toggle button").forEach(function (b) {
    b.addEventListener("click", function () { apply(b.dataset.themeSet); });
  });

  // follow OS changes while in "system" mode
  if (window.matchMedia) {
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", function () {
      if ((html.getAttribute("data-theme") || "system") === "system") apply("system");
    });
  }

  apply(current);
})();
