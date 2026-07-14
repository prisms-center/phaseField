function copy_parameter(element) {
    navigator.clipboard.writeText(element.textContent);

    element.classList.add("copied");

    setTimeout(() => {
        element.classList.remove("copied");
    }, 500);
}
