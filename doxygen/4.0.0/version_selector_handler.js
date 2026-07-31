$(function () {
    var repoName = window.location.pathname.split('/')[1];
    var doxygenFolderName = window.location.pathname.split('/')[2];

    var filePath = '/' + repoName + '/' + doxygenFolderName;

    $.get(filePath + '/version_selector.html', function (data) {
        // Inject version selector HTML into the page
        $('#projectnumber').html(data);

        // Event listener to handle version selection
        document.getElementById('versionSelector').addEventListener('change', function () {
            var selectedVersion = this.value;
            window.location.href = filePath + '/' + selectedVersion + '/index.html';
        });

        // Set the selected option based on the current version
        var currentVersion = window.location.pathname.split('/')[3];
        $('#versionSelector').val(currentVersion);
    });
});
