var host = window.location.host;
var isLocal = !host.includes('47.120.55.223')

var instructionsBtn = document.getElementById('instructions_title');
var flat = true
instructionsBtn.addEventListener('click', function () {
	console.log(flat)
	document.getElementById('instructions').style.opacity = flat ? 1 : 0
	flat = !flat
});