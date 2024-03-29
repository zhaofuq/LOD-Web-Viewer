
const settings = {
	renderingMode: "Octree", // "Octree", "Gaussian"
	renderResolution: 1.0,
	maxGaussians: 2000000,
	scalingBeta: 2,
	bgColor: '#000000',
	speed: 0.07,
	fov: 47,
	lodLevel: 4,
	en_frustum_culling: false,
	debugDepth: false,
	sortTime: 'NaN',
	fps: 0,
}
let activeKeys = [];

let cameras = [
	{
		id: 0,
		img_name: "00001",
		width: 1959,
		height: 1090,
		position: [
			-3.0089893469241797, -0.11086489695181866, -3.7527640949141428,
		],
		rotation: [
			[0.876134201218856, 0.06925962026449776, 0.47706599800804744],
			[-0.04747421839895102, 0.9972110940209488, -0.057586739349882114],
			[-0.4797239414934443, 0.027805376500959853, 0.8769787916452908],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 1,
		img_name: "00009",
		width: 1959,
		height: 1090,
		position: [
			-2.5199776022057296, -0.09704735754873686, -3.6247725540304545,
		],
		rotation: [
			[0.9982731285632193, -0.011928707708098955, -0.05751927260507243],
			[0.0065061360949636325, 0.9955928229282383, -0.09355533724430458],
			[0.058381769258182864, 0.09301955098900708, 0.9939511719154457],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 2,
		img_name: "00017",
		width: 1959,
		height: 1090,
		position: [
			-0.7737533667465242, -0.3364271945329695, -2.9358969417573753,
		],
		rotation: [
			[0.9998813418672372, 0.013742375651625236, -0.0069605529394208224],
			[-0.014268370388586709, 0.996512943252834, -0.08220929105659476],
			[0.00580653013657589, 0.08229885200307129, 0.9965907801935302],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 3,
		img_name: "00025",
		width: 1959,
		height: 1090,
		position: [
			1.2198221749590001, -0.2196687861401182, -2.3183162007028453,
		],
		rotation: [
			[0.9208648867765482, 0.0012010625395201253, 0.389880004297208],
			[-0.06298204172269357, 0.987319521752825, 0.14571693239364383],
			[-0.3847611242348369, -0.1587410451475895, 0.9092635249821667],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 4,
		img_name: "00033",
		width: 1959,
		height: 1090,
		position: [
			1.742387858893817, -0.13848225198886954, -2.0566370113193146,
		],
		rotation: [
			[0.24669889292141334, -0.08370189346592856, -0.9654706879349405],
			[0.11343747891376445, 0.9919082664242816, -0.05700815184573074],
			[0.9624300466054861, -0.09545671285663988, 0.2541976029815521],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 5,
		img_name: "00041",
		width: 1959,
		height: 1090,
		position: [
			3.6567309419223935, -0.16470990600750707, -1.3458085590422042,
		],
		rotation: [
			[0.2341293058324528, -0.02968330457755884, -0.9717522161434825],
			[0.10270823606832301, 0.99469554638321, -0.005638106875665722],
			[0.9667649592295676, -0.09848690996657204, 0.2359360976431732],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 6,
		img_name: "00049",
		width: 1959,
		height: 1090,
		position: [
			3.9013554243203497, -0.2597500978038105, -0.8106154188297828,
		],
		rotation: [
			[0.6717235545638952, -0.015718162115524837, -0.7406351366386528],
			[0.055627354673906296, 0.9980224478387622, 0.029270992841185218],
			[0.7387104058127439, -0.060861588786650656, 0.6712695459756353],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 7,
		img_name: "00057",
		width: 1959,
		height: 1090,
		position: [4.742994605467533, -0.05591660945412069, 0.9500365976084458],
		rotation: [
			[-0.17042655709210375, 0.01207080756938, -0.9852964448542146],
			[0.1165090336695526, 0.9931575292530063, -0.00798543433078162],
			[0.9784581921120181, -0.1161568667478904, -0.1706667764862097],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 8,
		img_name: "00065",
		width: 1959,
		height: 1090,
		position: [4.34676307626522, 0.08168160516967145, 1.0876221470355405],
		rotation: [
			[-0.003575447631888379, -0.044792503246552894, -0.9989899137764799],
			[0.10770152645126597, 0.9931680875192705, -0.04491693593046672],
			[0.9941768441149182, -0.10775333677534978, 0.0012732004866391048],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
	{
		id: 9,
		img_name: "00073",
		width: 1959,
		height: 1090,
		position: [3.264984351114202, 0.078974937336732, 1.0117200284114904],
		rotation: [
			[-0.026919994628162257, -0.1565891128261527, -0.9872968974090509],
			[0.08444552208239385, 0.983768234577625, -0.1583319754069128],
			[0.9960643893290491, -0.0876350978794554, -0.013259786205163005],
		],
		fy: 1164.6601287484507,
		fx: 1159.5880733038064,
	},
];

let camera = cameras[0];
let reloadLod = false;

const startReloadLod = () => {
	setTimeout(() => {
		if (!reloadLod) {
			reloadLod = true
		}
		setTimeout(() => {
			reloadLod = false;
		}, 1000);
	}, 10);

};

function getProjectionMatrix(fx, fy, width, height) {
	const znear = 0.2;
	const zfar = 200;
	return [
		[(2 * fx) / width, 0, 0, 0],
		[0, -(2 * fy) / height, 0, 0],
		[0, 0, zfar / (zfar - znear), 1],
		[0, 0, -(zfar * znear) / (zfar - znear), 0],
	].flat();
}

function getViewMatrix(camera) {
	const R = camera.rotation.flat();
	const t = camera.position;
	const camToWorld = [
		[R[0], R[1], R[2], 0],
		[R[3], R[4], R[5], 0],
		[R[6], R[7], R[8], 0],
		[
			-t[0] * R[0] - t[1] * R[3] - t[2] * R[6],
			-t[0] * R[1] - t[1] * R[4] - t[2] * R[7],
			-t[0] * R[2] - t[1] * R[5] - t[2] * R[8],
			1,
		],
	].flat();
	return camToWorld;
}

function multiply4(a, b) {
	return [
		b[0] * a[0] + b[1] * a[4] + b[2] * a[8] + b[3] * a[12],
		b[0] * a[1] + b[1] * a[5] + b[2] * a[9] + b[3] * a[13],
		b[0] * a[2] + b[1] * a[6] + b[2] * a[10] + b[3] * a[14],
		b[0] * a[3] + b[1] * a[7] + b[2] * a[11] + b[3] * a[15],
		b[4] * a[0] + b[5] * a[4] + b[6] * a[8] + b[7] * a[12],
		b[4] * a[1] + b[5] * a[5] + b[6] * a[9] + b[7] * a[13],
		b[4] * a[2] + b[5] * a[6] + b[6] * a[10] + b[7] * a[14],
		b[4] * a[3] + b[5] * a[7] + b[6] * a[11] + b[7] * a[15],
		b[8] * a[0] + b[9] * a[4] + b[10] * a[8] + b[11] * a[12],
		b[8] * a[1] + b[9] * a[5] + b[10] * a[9] + b[11] * a[13],
		b[8] * a[2] + b[9] * a[6] + b[10] * a[10] + b[11] * a[14],
		b[8] * a[3] + b[9] * a[7] + b[10] * a[11] + b[11] * a[15],
		b[12] * a[0] + b[13] * a[4] + b[14] * a[8] + b[15] * a[12],
		b[12] * a[1] + b[13] * a[5] + b[14] * a[9] + b[15] * a[13],
		b[12] * a[2] + b[13] * a[6] + b[14] * a[10] + b[15] * a[14],
		b[12] * a[3] + b[13] * a[7] + b[14] * a[11] + b[15] * a[15],
	];
}

function invert4(a) {
	let b00 = a[0] * a[5] - a[1] * a[4];
	let b01 = a[0] * a[6] - a[2] * a[4];
	let b02 = a[0] * a[7] - a[3] * a[4];
	let b03 = a[1] * a[6] - a[2] * a[5];
	let b04 = a[1] * a[7] - a[3] * a[5];
	let b05 = a[2] * a[7] - a[3] * a[6];
	let b06 = a[8] * a[13] - a[9] * a[12];
	let b07 = a[8] * a[14] - a[10] * a[12];
	let b08 = a[8] * a[15] - a[11] * a[12];
	let b09 = a[9] * a[14] - a[10] * a[13];
	let b10 = a[9] * a[15] - a[11] * a[13];
	let b11 = a[10] * a[15] - a[11] * a[14];
	let det =
		b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
	if (!det) return null;
	return [
		(a[5] * b11 - a[6] * b10 + a[7] * b09) / det,
		(a[2] * b10 - a[1] * b11 - a[3] * b09) / det,
		(a[13] * b05 - a[14] * b04 + a[15] * b03) / det,
		(a[10] * b04 - a[9] * b05 - a[11] * b03) / det,
		(a[6] * b08 - a[4] * b11 - a[7] * b07) / det,
		(a[0] * b11 - a[2] * b08 + a[3] * b07) / det,
		(a[14] * b02 - a[12] * b05 - a[15] * b01) / det,
		(a[8] * b05 - a[10] * b02 + a[11] * b01) / det,
		(a[4] * b10 - a[5] * b08 + a[7] * b06) / det,
		(a[1] * b08 - a[0] * b10 - a[3] * b06) / det,
		(a[12] * b04 - a[13] * b02 + a[15] * b00) / det,
		(a[9] * b02 - a[8] * b04 - a[11] * b00) / det,
		(a[5] * b07 - a[4] * b09 - a[6] * b06) / det,
		(a[0] * b09 - a[1] * b07 + a[2] * b06) / det,
		(a[13] * b01 - a[12] * b03 - a[14] * b00) / det,
		(a[8] * b03 - a[9] * b01 + a[10] * b00) / det,
	];
}

function rotate4(a, rad, x, y, z) {
	let len = Math.hypot(x, y, z);
	x /= len;
	y /= len;
	z /= len;
	let s = Math.sin(rad);
	let c = Math.cos(rad);
	let t = 1 - c;
	let b00 = x * x * t + c;
	let b01 = y * x * t + z * s;
	let b02 = z * x * t - y * s;
	let b10 = x * y * t - z * s;
	let b11 = y * y * t + c;
	let b12 = z * y * t + x * s;
	let b20 = x * z * t + y * s;
	let b21 = y * z * t - x * s;
	let b22 = z * z * t + c;
	return [
		a[0] * b00 + a[4] * b01 + a[8] * b02,
		a[1] * b00 + a[5] * b01 + a[9] * b02,
		a[2] * b00 + a[6] * b01 + a[10] * b02,
		a[3] * b00 + a[7] * b01 + a[11] * b02,
		a[0] * b10 + a[4] * b11 + a[8] * b12,
		a[1] * b10 + a[5] * b11 + a[9] * b12,
		a[2] * b10 + a[6] * b11 + a[10] * b12,
		a[3] * b10 + a[7] * b11 + a[11] * b12,
		a[0] * b20 + a[4] * b21 + a[8] * b22,
		a[1] * b20 + a[5] * b21 + a[9] * b22,
		a[2] * b20 + a[6] * b21 + a[10] * b22,
		a[3] * b20 + a[7] * b21 + a[11] * b22,
		...a.slice(12, 16),
	];
}

function translate4(a, x, y, z) {
	return [
		...a.slice(0, 12),
		a[0] * x + a[4] * y + a[8] * z + a[12],
		a[1] * x + a[5] * y + a[9] * z + a[13],
		a[2] * x + a[6] * y + a[10] * z + a[14],
		a[3] * x + a[7] * y + a[11] * z + a[15],
	];
}

function createWorker(self) {
	let buffer;
	let vertexCount = 0;
	let baseVertexCount = 0;

	let viewProj;
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	let lastProj = [];
	let depthIndex = new Uint32Array();
	let lastVertexCount = 0;

	var _floatView = new Float32Array(1);
	var _int32View = new Int32Array(_floatView.buffer);

	function floatToHalf(float) {
		_floatView[0] = float;
		var f = _int32View[0];

		var sign = (f >> 31) & 0x0001;
		var exp = (f >> 23) & 0x00ff;
		var frac = f & 0x007fffff;

		var newExp;
		if (exp == 0) {
			newExp = 0;
		} else if (exp < 113) {
			newExp = 0;
			frac |= 0x00800000;
			frac = frac >> (113 - exp);
			if (frac & 0x01000000) {
				newExp = 1;
				frac = 0;
			}
		} else if (exp < 142) {
			newExp = exp - 112;
		} else {
			newExp = 31;
			frac = 0;
		}

		return (sign << 15) | (newExp << 10) | (frac >> 13);
	}

	function packHalf2x16(x, y) {
		return (floatToHalf(x) | (floatToHalf(y) << 16)) >>> 0;
	}

	function generateTexture() {
		if (!buffer) return;
		const f_buffer = new Float32Array(buffer);
		const u_buffer = new Uint8Array(buffer);

		var texwidth = 1024 * 2; // Set to your desired width
		var texheight = Math.ceil((2 * vertexCount) / texwidth); // Set to your desired height
		var texdata = new Uint32Array(texwidth * texheight * 4); // 4 components per pixel (RGBA)
		var texdata_c = new Uint8Array(texdata.buffer);
		var texdata_f = new Float32Array(texdata.buffer);

		// Here we convert from a .splat file buffer into a texture
		// With a little bit more foresight perhaps this texture file
		// should have been the native format as it'd be very easy to
		// load it into webgl.
		for (let i = 0; i < vertexCount; i++) {
			// x, y, z
			texdata_f[8 * i + 0] = f_buffer[8 * i + 0];
			texdata_f[8 * i + 1] = f_buffer[8 * i + 1];
			texdata_f[8 * i + 2] = f_buffer[8 * i + 2];

			// r, g, b, a
			texdata_c[4 * (8 * i + 7) + 0] = u_buffer[32 * i + 24 + 0];
			texdata_c[4 * (8 * i + 7) + 1] = u_buffer[32 * i + 24 + 1];
			texdata_c[4 * (8 * i + 7) + 2] = u_buffer[32 * i + 24 + 2];
			texdata_c[4 * (8 * i + 7) + 3] = u_buffer[32 * i + 24 + 3];

			// quaternions
			let scale = [
				f_buffer[8 * i + 3 + 0],
				f_buffer[8 * i + 3 + 1],
				f_buffer[8 * i + 3 + 2],
			];
			let rot = [
				(u_buffer[32 * i + 28 + 0] - 128) / 128,
				(u_buffer[32 * i + 28 + 1] - 128) / 128,
				(u_buffer[32 * i + 28 + 2] - 128) / 128,
				(u_buffer[32 * i + 28 + 3] - 128) / 128,
			];

			// Compute the matrix product of S and R (M = S * R)
			const M = [
				1.0 - 2.0 * (rot[2] * rot[2] + rot[3] * rot[3]),
				2.0 * (rot[1] * rot[2] + rot[0] * rot[3]),
				2.0 * (rot[1] * rot[3] - rot[0] * rot[2]),

				2.0 * (rot[1] * rot[2] - rot[0] * rot[3]),
				1.0 - 2.0 * (rot[1] * rot[1] + rot[3] * rot[3]),
				2.0 * (rot[2] * rot[3] + rot[0] * rot[1]),

				2.0 * (rot[1] * rot[3] + rot[0] * rot[2]),
				2.0 * (rot[2] * rot[3] - rot[0] * rot[1]),
				1.0 - 2.0 * (rot[1] * rot[1] + rot[2] * rot[2]),
			].map((k, i) => k * scale[Math.floor(i / 3)]);

			const sigma = [
				M[0] * M[0] + M[3] * M[3] + M[6] * M[6],
				M[0] * M[1] + M[3] * M[4] + M[6] * M[7],
				M[0] * M[2] + M[3] * M[5] + M[6] * M[8],
				M[1] * M[1] + M[4] * M[4] + M[7] * M[7],
				M[1] * M[2] + M[4] * M[5] + M[7] * M[8],
				M[2] * M[2] + M[5] * M[5] + M[8] * M[8],
			];

			texdata[8 * i + 4] = packHalf2x16(4 * sigma[0], 4 * sigma[1]);
			texdata[8 * i + 5] = packHalf2x16(4 * sigma[2], 4 * sigma[3]);
			texdata[8 * i + 6] = packHalf2x16(4 * sigma[4], 4 * sigma[5]);
		}

		self.postMessage({ texdata, texwidth, texheight }, [texdata.buffer]);
	}

	function runSort(viewProj) {
		if (!buffer) return;
		const f_buffer = new Float32Array(buffer);
		if (lastVertexCount == vertexCount) {
			let dot =
				lastProj[2] * viewProj[2] +
				lastProj[6] * viewProj[6] +
				lastProj[10] * viewProj[10];
			if (Math.abs(dot - 1) < 0.01) {
				return;
			}
		} else {
			generateTexture();
			lastVertexCount = vertexCount;
		}

		const start = performance.now()
		let maxDepth = -Infinity;
		let minDepth = Infinity;
		let sizeList = new Int32Array(vertexCount);
		for (let i = 0; i < vertexCount; i++) {
			let depth =
				((viewProj[2] * f_buffer[8 * i + 0] +
					viewProj[6] * f_buffer[8 * i + 1] +
					viewProj[10] * f_buffer[8 * i + 2]) *
					4096) |
				0;
			sizeList[i] = depth;
			if (depth > maxDepth) maxDepth = depth;
			if (depth < minDepth) minDepth = depth;
		}

		// This is a 16 bit single-pass counting sort
		let depthInv = (256 * 256) / (maxDepth - minDepth);
		let counts0 = new Uint32Array(256 * 256);
		for (let i = 0; i < vertexCount; i++) {
			sizeList[i] = ((sizeList[i] - minDepth) * depthInv) | 0;
			counts0[sizeList[i]]++;
		}
		let starts0 = new Uint32Array(256 * 256);
		for (let i = 1; i < 256 * 256; i++)
			starts0[i] = starts0[i - 1] + counts0[i - 1];
		depthIndex = new Uint32Array(vertexCount);
		for (let i = 0; i < vertexCount; i++)
			depthIndex[starts0[sizeList[i]]++] = i;

		const sortTime = `${((performance.now() - start) / 1000).toFixed(3)}s`

		lastProj = viewProj;
		self.postMessage({ depthIndex, viewProj, vertexCount, sortTime }, [
			depthIndex.buffer,
		]);
	}

	function processPlyBuffer(inputBuffer) {
		const ubuf = new Uint8Array(inputBuffer);
		// 10KB ought to be enough for a header...
		const header = new TextDecoder().decode(ubuf.slice(0, 1024 * 10));
		const header_end = "end_header\n";
		const header_end_index = header.indexOf(header_end);
		if (header_end_index < 0)
			throw new Error("Unable to read .ply file header");
		const vertexCount = parseInt(/element vertex (\d+)\n/.exec(header)[1]);
		console.log("Vertex Count", vertexCount);
		let row_offset = 0,
			offsets = {},
			types = {};
		const TYPE_MAP = {
			double: "getFloat64",
			int: "getInt32",
			uint: "getUint32",
			float: "getFloat32",
			short: "getInt16",
			ushort: "getUint16",
			uchar: "getUint8",
		};
		for (let prop of header
			.slice(0, header_end_index)
			.split("\n")
			.filter((k) => k.startsWith("property "))) {
			const [p, type, name] = prop.split(" ");
			const arrayType = TYPE_MAP[type] || "getInt8";
			types[name] = arrayType;
			offsets[name] = row_offset;
			row_offset += parseInt(arrayType.replace(/[^\d]/g, "")) / 8;
		}
		console.log("Bytes per row", row_offset, types, offsets);

		let dataView = new DataView(
			inputBuffer,
			header_end_index + header_end.length,
		);
		let row = 0;
		const attrs = new Proxy(
			{},
			{
				get(target, prop) {
					if (!types[prop]) throw new Error(prop + " not found");
					return dataView[types[prop]](
						row * row_offset + offsets[prop],
						true,
					);
				},
			},
		);

		console.time("calculate importance");
		let sizeList = new Float32Array(vertexCount);
		let sizeIndex = new Uint32Array(vertexCount);
		for (row = 0; row < vertexCount; row++) {
			sizeIndex[row] = row;
			if (!types["scale_0"]) continue;
			const size =
				Math.exp(attrs.scale_0) *
				Math.exp(attrs.scale_1) *
				Math.exp(attrs.scale_2);
			const opacity = 1 / (1 + Math.exp(-attrs.opacity));
			sizeList[row] = size * opacity;
		}
		console.timeEnd("calculate importance");

		console.time("sort");
		sizeIndex.sort((b, a) => sizeList[a] - sizeList[b]);
		console.timeEnd("sort");

		// 6*4 + 4 + 4 = 8*4
		// XYZ - Position (Float32)
		// XYZ - Scale (Float32)
		// RGBA - colors (uint8)
		// IJKL - quaternion/rot (uint8)
		const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
		const buffer = new ArrayBuffer(rowLength * vertexCount);

		console.time("build buffer");
		for (let j = 0; j < vertexCount; j++) {
			row = sizeIndex[j];

			const position = new Float32Array(buffer, j * rowLength, 3);
			const scales = new Float32Array(buffer, j * rowLength + 4 * 3, 3);
			const rgba = new Uint8ClampedArray(
				buffer,
				j * rowLength + 4 * 3 + 4 * 3,
				4,
			);
			const rot = new Uint8ClampedArray(
				buffer,
				j * rowLength + 4 * 3 + 4 * 3 + 4,
				4,
			);

			if (types["scale_0"]) {
				const qlen = Math.sqrt(
					attrs.rot_0 ** 2 +
					attrs.rot_1 ** 2 +
					attrs.rot_2 ** 2 +
					attrs.rot_3 ** 2,
				);

				rot[0] = (attrs.rot_0 / qlen) * 128 + 128;
				rot[1] = (attrs.rot_1 / qlen) * 128 + 128;
				rot[2] = (attrs.rot_2 / qlen) * 128 + 128;
				rot[3] = (attrs.rot_3 / qlen) * 128 + 128;

				scales[0] = Math.exp(attrs.scale_0);
				scales[1] = Math.exp(attrs.scale_1);
				scales[2] = Math.exp(attrs.scale_2);
			} else {
				scales[0] = 0.01;
				scales[1] = 0.01;
				scales[2] = 0.01;

				rot[0] = 255;
				rot[1] = 0;
				rot[2] = 0;
				rot[3] = 0;
			}

			position[0] = attrs.x;
			position[1] = attrs.y;
			position[2] = attrs.z;

			if (types["f_dc_0"]) {
				const SH_C0 = 0.28209479177387814;
				rgba[0] = (0.5 + SH_C0 * attrs.f_dc_0) * 255;
				rgba[1] = (0.5 + SH_C0 * attrs.f_dc_1) * 255;
				rgba[2] = (0.5 + SH_C0 * attrs.f_dc_2) * 255;
			} else {
				rgba[0] = attrs.red;
				rgba[1] = attrs.green;
				rgba[2] = attrs.blue;
			}
			if (types["opacity"]) {
				rgba[3] = (1 / (1 + Math.exp(-attrs.opacity))) * 255;
			} else {
				rgba[3] = 255;
			}
		}
		console.timeEnd("build buffer");
		return buffer;
	}

	const throttledSort = () => {
		if (!sortRunning) {
			sortRunning = true;
			let lastView = viewProj;
			runSort(lastView);
			setTimeout(() => {
				sortRunning = false;
				if (lastView !== viewProj) {
					throttledSort();
				}
			}, 0);
		}
	};

	let sortRunning;
	self.onmessage = (e) => {
		if (e.data.ply) {
			vertexCount = 0;
			runSort(viewProj);
			buffer = processPlyBuffer(e.data.ply);
			vertexCount = Math.floor(buffer.byteLength / rowLength);
			postMessage({ buffer: buffer });
		} else if (e.data.buffer) {
			buffer = e.data.buffer;
			vertexCount = e.data.vertexCount;
			console.log("Recive Vertex Count", vertexCount);
		} else if (e.data.vertexCount) {
			vertexCount = e.data.vertexCount;
		} else if (e.data.view) {
			viewProj = e.data.view;
			throttledSort();
		}
	};
}

const vertexShaderSource = `
#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;

in vec2 position;
in int index;

out vec4 vColor;
out vec2 vPosition;

void main () {
    uvec4 cen = texelFetch(u_texture, ivec2((uint(index) & 0x3ffu) << 1, uint(index) >> 10), 0);
    vec4 cam = view * vec4(uintBitsToFloat(cen.xyz), 1);
    vec4 pos2d = projection * cam;

    float clip = 1.2 * pos2d.w;
    if (pos2d.z < -clip || pos2d.x < -clip || pos2d.x > clip || pos2d.y < -clip || pos2d.y > clip) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

    uvec4 cov = texelFetch(u_texture, ivec2(((uint(index) & 0x3ffu) << 1) | 1u, uint(index) >> 10), 0);
    vec2 u1 = unpackHalf2x16(cov.x), u2 = unpackHalf2x16(cov.y), u3 = unpackHalf2x16(cov.z);
    mat3 Vrk = mat3(u1.x, u1.y, u2.x, u1.y, u2.y, u3.x, u2.x, u3.x, u3.y);

    mat3 J = mat3(
        focal.x / cam.z, 0., -(focal.x * cam.x) / (cam.z * cam.z), 
        0., -focal.y / cam.z, (focal.y * cam.y) / (cam.z * cam.z), 
        0., 0., 0.
    );

    mat3 T = transpose(mat3(view)) * J;
    mat3 cov2d = transpose(T) * Vrk * T;

    float mid = (cov2d[0][0] + cov2d[1][1]) / 2.0;
    float radius = length(vec2((cov2d[0][0] - cov2d[1][1]) / 2.0, cov2d[0][1]));
    float lambda1 = mid + radius, lambda2 = mid - radius;

    if(lambda2 < 0.0) return;
    vec2 diagonalVector = normalize(vec2(cov2d[0][1], lambda1 - cov2d[0][0]));
    vec2 majorAxis = min(sqrt(2.0 * lambda1), 1024.0) * diagonalVector;
    vec2 minorAxis = min(sqrt(2.0 * lambda2), 1024.0) * vec2(diagonalVector.y, -diagonalVector.x);

    vColor = clamp(pos2d.z/pos2d.w+1.0, 0.0, 1.0) * vec4((cov.w) & 0xffu, (cov.w >> 8) & 0xffu, (cov.w >> 16) & 0xffu, (cov.w >> 24) & 0xffu) / 255.0;
    vPosition = position;

    vec2 vCenter = vec2(pos2d) / pos2d.w;
    gl_Position = vec4(
        vCenter 
        + position.x * majorAxis / viewport 
        + position.y * minorAxis / viewport, 0.0, 1.0);

}
`.trim();

const fragmentShaderSource = `
#version 300 es
precision highp float;

in vec4 vColor;
in vec2 vPosition;

out vec4 fragColor;

void main () {
    float A = -dot(vPosition, vPosition);
    if (A < -4.0) discard;
    float B = exp(A) * vColor.a;
    fragColor = vec4(B * vColor.rgb, B);
}

`.trim();

// let defaultViewMatrix = [0.578, 0.190, -0.790, 0, 
//                         0.814, -0.142, 0.555, 0, 
//                         -0.011, -0.975, -0.237, 0, 
//                         0.104, 0.080, 3.48, 0.999]
let defaultViewMatrix = dataSource.defaultViewMatrix

let viewMatrix = defaultViewMatrix;
let lastViewMatrix = viewMatrix;
let projectionMatrix;
let stopReading = false;
let worker = null;
let octreeGeometry, octreeGeometryLoader;
const gaussianSplats = {}
import { loadPointCloud } from "./potree.js"

// Init settings GUI panel
function initGUI(resize) {
	const gui = new lil.GUI({ title: 'LetsGo Settings' })
	gui.add(settings, 'renderingMode', ['Gaussian', 'Octree']).name('Rendering Mode')
	gui.add(settings, 'sortTime').name('Sort Time').disable().listen()
	gui.add(settings, 'fps').name('FPS').disable().listen()
	// gui.add(settings, 'renderResolution', 0.1, 1, 0.01).name('Preview Resolution')
	//     .onChange(() =>  resize())



	gui.add(settings, 'lodLevel', 4, 16, 1).name('LOD Level')
		.onChange(async value => {
			try {
				stopReading = true;
				await updateGaussianByView(viewMatrix, projectionMatrix, value, settings.maxGaussians)
			} catch (error) {
				throw error
			}
		})

	gui.add(settings, 'maxGaussians', 1e4, 1e7, 10000).name('Ext Gaussians')
		.onChange(async value => {
			try {
				stopReading = true;
				await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, value)
			} catch (error) {
				throw error
			}
		})

	gui.add(settings, 'scalingBeta', 1, 16, 1).name('scalingBeta')
		.onChange(async value => {
			try {
				stopReading = true;
				await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, value)
			} catch (error) {
				throw error
			}
		})


	// 创建一个包含说明文案的HTML元素
	const description = document.createElement('div');
	description.innerHTML = ' Parameters:\
	<br /> LOD Level: Controls the level-of-detail of 3DGS model. Higher values result in more details \
	<br /><br /> Ext Gaussians: Controls the extra number of Gaussians when level > 4. Increasing this value adds detail and complexity to the scene, but may also increase loading times and performance overhead. Adjust based on scene complexity and device capabilities. \
	<br /><br /> scalingBeta: Controls coef of the level decreasing from distance. Increasing this value will result in fast level decreasing from neaf to far';
	description.style.fontSize = '12px';
	description.style.color = '#fff';
	description.style.padding = '8px';
	description.style.opacity = 0.8;
	description.style.lineHeight = '150%';

	// 将HTML元素添加到GUI界面中
	const firstChildDiv = gui.domElement.querySelector('.children');
	firstChildDiv.appendChild(description);


}
const chunksData = [];
async function main() {
	let carousel = true;

	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	let splatData = new Uint8Array(rowLength * 1000000);
	// const downsample =
	//     splatData.length / rowLength > 500000 ? 1 : 1 / devicePixelRatio;
	const downsample = settings.renderResolution;
	// Worker Init settings
	worker = new Worker(
		URL.createObjectURL(
			new Blob(["(", createWorker.toString(), ")(self)"], {
				type: "application/javascript",
			}),
		),
	);

	// GL setup
	const canvas = document.getElementById("canvas");

	const gl = canvas.getContext("webgl2", {
		antialias: false,
	});

	const vertexShader = gl.createShader(gl.VERTEX_SHADER);
	gl.shaderSource(vertexShader, vertexShaderSource);
	gl.compileShader(vertexShader);
	if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS))
		console.error(gl.getShaderInfoLog(vertexShader));

	const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
	gl.shaderSource(fragmentShader, fragmentShaderSource);
	gl.compileShader(fragmentShader);
	if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS))
		console.error(gl.getShaderInfoLog(fragmentShader));

	const program = gl.createProgram();
	gl.attachShader(program, vertexShader);
	gl.attachShader(program, fragmentShader);
	gl.linkProgram(program);
	gl.useProgram(program);

	if (!gl.getProgramParameter(program, gl.LINK_STATUS))
		console.error(gl.getProgramInfoLog(program));

	gl.disable(gl.DEPTH_TEST); // Disable depth testing

	// Enable blending
	gl.enable(gl.BLEND);
	gl.blendFuncSeparate(
		gl.ONE_MINUS_DST_ALPHA,
		gl.ONE,
		gl.ONE_MINUS_DST_ALPHA,
		gl.ONE,
	);
	gl.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);

	const u_projection = gl.getUniformLocation(program, "projection");
	const u_viewport = gl.getUniformLocation(program, "viewport");
	const u_focal = gl.getUniformLocation(program, "focal");
	const u_view = gl.getUniformLocation(program, "view");

	// positions
	const triangleVertices = new Float32Array([-2, -2, 2, -2, 2, 2, -2, 2]);
	const vertexBuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);
	const a_position = gl.getAttribLocation(program, "position");
	gl.enableVertexAttribArray(a_position);
	gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
	gl.vertexAttribPointer(a_position, 2, gl.FLOAT, false, 0, 0);

	var texture = gl.createTexture();
	gl.bindTexture(gl.TEXTURE_2D, texture);

	var u_textureLocation = gl.getUniformLocation(program, "u_texture");
	gl.uniform1i(u_textureLocation, 0);

	const indexBuffer = gl.createBuffer();
	const a_index = gl.getAttribLocation(program, "index");
	gl.enableVertexAttribArray(a_index);
	gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
	gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
	gl.vertexAttribDivisor(a_index, 1);

	const resize = () => {
		gl.uniform2fv(u_focal, new Float32Array([camera.fx, camera.fy]));

		projectionMatrix = getProjectionMatrix(
			camera.fx,
			camera.fy,
			innerWidth,
			innerHeight,
		);

		gl.uniform2fv(u_viewport, new Float32Array([innerWidth, innerHeight]));

		gl.canvas.width = Math.round(innerWidth / downsample);
		gl.canvas.height = Math.round(innerHeight / downsample);
		gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

		gl.uniformMatrix4fv(u_projection, false, projectionMatrix);
	};

	window.addEventListener("resize", resize);
	resize();

	worker.onmessage = (e) => {
		if (e.data.buffer) {
			splatData = new Uint8Array(e.data.buffer);
			const blob = new Blob([splatData.buffer], {
				type: "application/octet-stream",
			});
			const link = document.createElement("a");
			link.download = "model.splat";
			link.href = URL.createObjectURL(blob);
			document.body.appendChild(link);
			link.click();
		} else if (e.data.texdata) {
			const { texdata, texwidth, texheight } = e.data;
			// console.log(texdata)
			gl.bindTexture(gl.TEXTURE_2D, texture);
			gl.texParameteri(
				gl.TEXTURE_2D,
				gl.TEXTURE_WRAP_S,
				gl.CLAMP_TO_EDGE,
			);
			gl.texParameteri(
				gl.TEXTURE_2D,
				gl.TEXTURE_WRAP_T,
				gl.CLAMP_TO_EDGE,
			);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

			gl.texImage2D(
				gl.TEXTURE_2D,
				0,
				gl.RGBA32UI,
				texwidth,
				texheight,
				0,
				gl.RGBA_INTEGER,
				gl.UNSIGNED_INT,
				texdata,
			);
			gl.activeTexture(gl.TEXTURE0);
			gl.bindTexture(gl.TEXTURE_2D, texture);
		} else if (e.data.depthIndex) {
			const { depthIndex, viewProj } = e.data;
			gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);
			vertexCount = e.data.vertexCount;
			settings.sortTime = e.data.sortTime;
		}
	};



	window.addEventListener("keydown", (e) => {
		// if (document.activeElement != document.body) return;
		carousel = false;
		if (!activeKeys.includes(e.code)) activeKeys.push(e.code);
		if (/\d/.test(e.key)) {
			camera = cameras[parseInt(e.key)];
			viewMatrix = getViewMatrix(camera);
		}
		if (e.code == "KeyV") {
			location.hash =
				"#" +
				JSON.stringify(
					viewMatrix.map((k) => Math.round(k * 100) / 100),
				);
		} else if (e.code === "KeyP") {
			carousel = true;
		}
	});
	window.addEventListener("keyup", (e) => {
		activeKeys = activeKeys.filter((k) => k !== e.code);
	});
	window.addEventListener("blur", () => {
		activeKeys = [];
	});

	window.addEventListener(
		"wheel",
		(e) => {
			carousel = false;
			e.preventDefault();
			const lineHeight = 10;
			const scale =
				e.deltaMode == 1
					? lineHeight
					: e.deltaMode == 2
						? innerHeight
						: 1;
			let inv = invert4(viewMatrix);
			if (e.shiftKey) {
				inv = translate4(
					inv,
					(e.deltaX * scale) / innerWidth,
					(e.deltaY * scale) / innerHeight,
					0,
				);
			} else if (e.ctrlKey || e.metaKey) {
				// inv = rotate4(inv,  (e.deltaX * scale) / innerWidth,  0, 0, 1);
				// inv = translate4(inv,  0, (e.deltaY * scale) / innerHeight, 0);
				// let preY = inv[13];
				inv = translate4(
					inv,
					0,
					0,
					(-10 * (e.deltaY * scale)) / innerHeight,
				);
				// inv[13] = preY;
			} else {
				let d = 4;
				inv = translate4(inv, 0, 0, d);
				inv = rotate4(inv, -(e.deltaX * scale) / innerWidth, 0, 1, 0);
				inv = rotate4(inv, (e.deltaY * scale) / innerHeight, 1, 0, 0);
				inv = translate4(inv, 0, 0, -d);
			}

			viewMatrix = invert4(inv);
		},
		{ passive: false },
	);

	let startX, startY, down;
	canvas.addEventListener("mousedown", (e) => {
		carousel = false;
		e.preventDefault();
		startX = e.clientX;
		startY = e.clientY;
		down = e.ctrlKey || e.metaKey ? 2 : 1;
	});

	canvas.addEventListener("contextmenu", (e) => {
		carousel = false;
		e.preventDefault();
		startX = e.clientX;
		startY = e.clientY;
		down = 2;
	});

	canvas.addEventListener("mousemove", (e) => {
		e.preventDefault();
		if (down == 1) {
			let inv = invert4(viewMatrix);
			let dx = (5 * (e.clientX - startX)) / innerWidth;
			let dy = (5 * (e.clientY - startY)) / innerHeight;
			let d = 4;

			inv = translate4(inv, 0, 0, d);
			inv = rotate4(inv, dx, 0, 1, 0);
			inv = rotate4(inv, -dy, 1, 0, 0);
			inv = translate4(inv, 0, 0, -d);
			// let postAngle = Math.atan2(inv[0], inv[10])
			// inv = rotate4(inv, postAngle - preAngle, 0, 0, 1)
			// console.log(postAngle)
			viewMatrix = invert4(inv);

			startX = e.clientX;
			startY = e.clientY;
		} else if (down == 2) {
			let inv = invert4(viewMatrix);
			// inv = rotateY(inv, );
			// let preY = inv[13];
			inv = translate4(
				inv,
				(-10 * (e.clientX - startX)) / innerWidth,
				0,
				(10 * (e.clientY - startY)) / innerHeight,
			);
			// inv[13] = preY;
			viewMatrix = invert4(inv);

			startX = e.clientX;
			startY = e.clientY;
		}
	});
	canvas.addEventListener("mouseup", (e) => {
		e.preventDefault();
		startReloadLod()
		console.log(viewMatrix)
		down = false;
		startX = 0;
		startY = 0;

	});

	let altX = 0,
		altY = 0;
	canvas.addEventListener(
		"touchstart",
		(e) => {
			e.preventDefault();
			if (e.touches.length === 1) {
				carousel = false;
				startX = e.touches[0].clientX;
				startY = e.touches[0].clientY;
				down = 1;
			} else if (e.touches.length === 2) {
				// console.log('beep')
				carousel = false;
				startX = e.touches[0].clientX;
				altX = e.touches[1].clientX;
				startY = e.touches[0].clientY;
				altY = e.touches[1].clientY;
				down = 1;
			}
		},
		{ passive: false },
	);
	canvas.addEventListener(
		"touchmove",
		(e) => {
			e.preventDefault();
			if (e.touches.length === 1 && down) {
				let inv = invert4(viewMatrix);
				let dx = (4 * (e.touches[0].clientX - startX)) / innerWidth;
				let dy = (4 * (e.touches[0].clientY - startY)) / innerHeight;

				let d = 4;
				inv = translate4(inv, 0, 0, d);
				// inv = translate4(inv,  -x, -y, -z);
				// inv = translate4(inv,  x, y, z);
				inv = rotate4(inv, dx, 0, 1, 0);
				inv = rotate4(inv, -dy, 1, 0, 0);
				inv = translate4(inv, 0, 0, -d);

				viewMatrix = invert4(inv);

				startX = e.touches[0].clientX;
				startY = e.touches[0].clientY;
			} else if (e.touches.length === 2) {
				// alert('beep')
				const dtheta =
					Math.atan2(startY - altY, startX - altX) -
					Math.atan2(
						e.touches[0].clientY - e.touches[1].clientY,
						e.touches[0].clientX - e.touches[1].clientX,
					);
				const dscale =
					Math.hypot(startX - altX, startY - altY) /
					Math.hypot(
						e.touches[0].clientX - e.touches[1].clientX,
						e.touches[0].clientY - e.touches[1].clientY,
					);
				const dx =
					(e.touches[0].clientX +
						e.touches[1].clientX -
						(startX + altX)) /
					2;
				const dy =
					(e.touches[0].clientY +
						e.touches[1].clientY -
						(startY + altY)) /
					2;
				let inv = invert4(viewMatrix);
				// inv = translate4(inv,  0, 0, d);
				inv = rotate4(inv, dtheta, 0, 0, 1);

				inv = translate4(inv, -dx / innerWidth, -dy / innerHeight, 0);

				// let preY = inv[13];
				inv = translate4(inv, 0, 0, 3 * (1 - dscale));
				// inv[13] = preY;

				viewMatrix = invert4(inv);

				startX = e.touches[0].clientX;
				altX = e.touches[1].clientX;
				startY = e.touches[0].clientY;
				altY = e.touches[1].clientY;
			}
		},
		{ passive: false },
	);
	canvas.addEventListener(
		"touchend",
		(e) => {
			e.preventDefault();
			startReloadLod()
			down = false;
			startX = 0;
			startY = 0;
		},
		{ passive: false },
	);

	let jumpDelta = 0;
	let vertexCount = 0;

	let lastFrame = 0;
	let avgFps = 0;
	let start = 0;

	window.addEventListener("gamepadconnected", (e) => {
		const gp = navigator.getGamepads()[e.gamepad.index];
		console.log(
			`Gamepad connected at index ${gp.index}: ${gp.id}. It has ${gp.buttons.length} buttons and ${gp.axes.length} axes.`,
		);
	});
	window.addEventListener("gamepaddisconnected", (e) => {
		console.log("Gamepad disconnected");
	});

	let leftGamepadTrigger, rightGamepadTrigger;

	const frame = (now) => {
		let inv = invert4(viewMatrix);
		let shiftKey = activeKeys.includes("Shift") || activeKeys.includes("ShiftLeft") || activeKeys.includes("ShiftRight")

		if (activeKeys.includes("KeyW")) {
			if (shiftKey) {
				inv = translate4(inv, 0, -0.03, 0);
			} else {
				inv = translate4(inv, 0, 0, 0.1);
			}
		}
		if (activeKeys.includes("KeyS")) {
			if (shiftKey) {
				inv = translate4(inv, 0, 0.03, 0);
			} else {
				inv = translate4(inv, 0, 0, -0.1);
			}
		}
		if (activeKeys.includes("KeyA"))
			inv = translate4(inv, -0.03, 0, 0);
		//
		if (activeKeys.includes("KeyD"))
			inv = translate4(inv, 0.03, 0, 0);
		// inv = rotate4(inv, 0.01, 0, 1, 0);
		// if (activeKeys.includes("KeyW")) inv = rotate4(inv, 0.005, 1, 0, 0);
		// if (activeKeys.includes("KeyS")) inv = rotate4(inv, -0.005, 1, 0, 0);
		// if (activeKeys.includes("KeyA")) inv = rotate4(inv, -0.01, 0, 1, 0);
		// if (activeKeys.includes("KeyD")) inv = rotate4(inv, 0.01, 0, 1, 0);

		if (activeKeys.includes("KeyQ")) inv = rotate4(inv, 0.01, 0, 0, 1);
		if (activeKeys.includes("KeyE")) inv = rotate4(inv, -0.01, 0, 0, 1);

		if (activeKeys.includes("Space")) {
			inv = translate4(inv, 0, -0.05, 0);
			inv = rotate4(inv, -0.005, 1, 0, 0);
		}

		if (activeKeys.includes("ShiftLeft")) {
			inv = translate4(inv, 0, 0.05, 0);
			inv = rotate4(inv, 0.005, 1, 0, 0);
		}

		const gamepads = navigator.getGamepads ? navigator.getGamepads() : [];
		let isJumping = activeKeys.includes("Ctrl");
		for (let gamepad of gamepads) {
			if (!gamepad) continue;

			const axisThreshold = 0.1; // Threshold to detect when the axis is intentionally moved
			const moveSpeed = 0.06;
			const rotateSpeed = 0.02;

			// Assuming the left stick controls translation (axes 0 and 1)
			if (Math.abs(gamepad.axes[0]) > axisThreshold) {
				inv = translate4(inv, moveSpeed * gamepad.axes[0], 0, 0);
				carousel = false;
			}
			if (Math.abs(gamepad.axes[1]) > axisThreshold) {
				inv = translate4(inv, 0, 0, -moveSpeed * gamepad.axes[1]);
				carousel = false;
			}
			if (gamepad.buttons[12].pressed || gamepad.buttons[13].pressed) {
				inv = translate4(inv, 0, -moveSpeed * (gamepad.buttons[12].pressed - gamepad.buttons[13].pressed), 0);
				carousel = false;
			}

			if (gamepad.buttons[14].pressed || gamepad.buttons[15].pressed) {
				inv = translate4(inv, -moveSpeed * (gamepad.buttons[14].pressed - gamepad.buttons[15].pressed), 0, 0);
				carousel = false;
			}

			// Assuming the right stick controls rotation (axes 2 and 3)
			if (Math.abs(gamepad.axes[2]) > axisThreshold) {
				inv = rotate4(inv, rotateSpeed * gamepad.axes[2], 0, 1, 0);
				carousel = false;
			}
			if (Math.abs(gamepad.axes[3]) > axisThreshold) {
				inv = rotate4(inv, -rotateSpeed * gamepad.axes[3], 1, 0, 0);
				carousel = false;
			}

			let tiltAxis = gamepad.buttons[6].value - gamepad.buttons[7].value;
			if (Math.abs(tiltAxis) > axisThreshold) {
				inv = rotate4(inv, rotateSpeed * tiltAxis, 0, 0, 1);
				carousel = false;
			}
			if (gamepad.buttons[4].pressed && !leftGamepadTrigger) {
				camera = cameras[(cameras.indexOf(camera) + 1) % cameras.length]
				inv = invert4(getViewMatrix(camera));
				carousel = false;
			}
			if (gamepad.buttons[5].pressed && !rightGamepadTrigger) {
				camera = cameras[(cameras.indexOf(camera) + cameras.length - 1) % cameras.length]
				inv = invert4(getViewMatrix(camera));
				carousel = false;
			}
			leftGamepadTrigger = gamepad.buttons[4].pressed;
			rightGamepadTrigger = gamepad.buttons[5].pressed;
			if (gamepad.buttons[0].pressed) {
				isJumping = true;
				carousel = false;
			}
			if (gamepad.buttons[3].pressed) {
				carousel = true;
			}
		}

		if (
			["KeyJ", "KeyK", "KeyL", "KeyI"].some((k) => activeKeys.includes(k))
		) {
			let d = 4;
			inv = translate4(inv, 0, 0, d);
			inv = rotate4(
				inv,
				activeKeys.includes("KeyJ")
					? -0.05
					: activeKeys.includes("KeyL")
						? 0.05
						: 0,
				0,
				1,
				0,
			);
			inv = rotate4(
				inv,
				activeKeys.includes("KeyI")
					? 0.05
					: activeKeys.includes("KeyK")
						? -0.05
						: 0,
				1,
				0,
				0,
			);
			inv = translate4(inv, 0, 0, -d);
		}

		viewMatrix = invert4(inv);

		if (carousel) {
			let inv = invert4(defaultViewMatrix);

			const t = Math.sin((Date.now() - start) / 5000);
			inv = translate4(inv, 2.5 * t, 0, 6 * (1 - Math.cos(t)));
			inv = rotate4(inv, -0.6 * t, 0, 1, 0);

			viewMatrix = invert4(inv);
		}

		const viewProj = multiply4(projectionMatrix, viewMatrix);
		worker.postMessage({ view: viewProj });

		const lastViewProj = multiply4(projectionMatrix, lastViewMatrix);
		const viewProjNorm = Math.sqrt(viewProj[2] * viewProj[2] + viewProj[6] * viewProj[6] + viewProj[10] * viewProj[10]);
		const lastViewProjNorm = Math.sqrt(lastViewProj[2] * lastViewProj[2] + lastViewProj[6] * lastViewProj[6] + lastViewProj[10] * lastViewProj[10]);

		let dot =
			lastViewProj[2] / lastViewProjNorm * viewProj[2] / viewProjNorm +
			lastViewProj[6] / lastViewProjNorm * viewProj[6] / viewProjNorm +
			lastViewProj[10] / lastViewProjNorm * viewProj[10] / viewProjNorm;

		if ((reloadLod || activeKeys.includes("ControlLeft")) && settings.renderingMode == "Octree" && Math.abs(dot - 1) > 0.01) {
			console.log(viewMatrix)
			stopReading = true;
			updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
			lastViewMatrix = viewMatrix;
		}

		const currentFps = 1000 / (now - lastFrame) || 0;
		avgFps = avgFps * 0.9 + currentFps * 0.1;

		if (vertexCount > 0) {
			// document.getElementById("spinner").style.display = "none";
			gl.uniformMatrix4fv(u_view, false, viewMatrix);
			gl.clear(gl.COLOR_BUFFER_BIT);
			gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);
		} else {
			gl.clear(gl.COLOR_BUFFER_BIT);
			// document.getElementById("spinner").style.display = "";
			start = Date.now() + 2000;
		}
		const progress = (100 * vertexCount) / (splatData.length / rowLength);
		if (progress < 95) {
			document.getElementById("progress").style.width = progress + "%";
		} else {
			document.getElementById("progress").style.display = "none";
			document.getElementById("spinner").style.display = "none";
		}

		settings.fps = Math.round(avgFps)
		lastFrame = now;
		requestAnimationFrame(frame);
	};

	frame();

	const selectFile = (file) => {
		const fr = new FileReader();

		// change render mode to gaussian-lod
		settings.mode = "gaussian";

		if (/\.json$/i.test(file.name)) {
			fr.onload = () => {
				cameras = JSON.parse(fr.result);
				viewMatrix = getViewMatrix(cameras[0]);
				projectionMatrix = getProjectionMatrix(
					camera.fx / downsample,
					camera.fy / downsample,
					canvas.width,
					canvas.height,
				);
				gl.uniformMatrix4fv(u_projection, false, projectionMatrix);

				console.log("Loaded Cameras");
			};
			fr.readAsText(file);
		} else {
			stopLoading = true;
			fr.onload = () => {
				splatData = new Uint8Array(fr.result);
				console.log("Loaded", Math.floor(splatData.length / rowLength));

				if (
					splatData[0] == 112 &&
					splatData[1] == 108 &&
					splatData[2] == 121 &&
					splatData[3] == 10
				) {
					// ply file magic header means it should be handled differently
					worker.postMessage({ ply: splatData.buffer });
				} else {
					worker.postMessage({
						buffer: splatData.buffer,
						vertexCount: Math.floor(splatData.length / rowLength),
					});
				}
			};
			fr.readAsArrayBuffer(file);
		}
	};

	window.addEventListener("hashchange", (e) => {
		try {
			viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
			carousel = false;
		} catch (err) { }
	});

	const preventDefault = (e) => {
		e.preventDefault();
		e.stopPropagation();
	};
	document.addEventListener("dragenter", preventDefault);
	document.addEventListener("dragover", preventDefault);

	// Setup Setting GUI
	initGUI(resize)

	// Init Gaussian Scene
	let bytesRead = 0;
	let lastVertexCount = -1;
	let stopLoading = false;

	// while (true) {
	// 	const { done, value } = await reader.read();
	// 	if (done || stopLoading) break;

	// 	splatData.set(value, bytesRead);
	// 	bytesRead += value.length;

	// 	if (vertexCount > lastVertexCount) {
	// 		// worker.postMessage({
	// 		// 	buffer: splatData.buffer,
	// 		// 	vertexCount: Math.floor(bytesRead / rowLength),
	// 		// });
	// 		lastVertexCount = vertexCount;
	// 	}
	// }
	// if (!stopLoading)
	// 	worker.postMessage({
	// 		buffer: splatData.buffer,
	// 		vertexCount: Math.floor(bytesRead / rowLength),
	// 	});
}

main().catch((err) => {
	// document.getElementById("spinner").style.display = "none";
	document.getElementById("message").innerText = err.toString();
});


function nodeView(i, node) {
	const att = node.geometry.attributes
	const min = node.boundingBox.min
	let position = [att.position.array[i * 3] + min.x,
	att.position.array[i * 3 + 1] + min.y,
	att.position.array[i * 3 + 2] + min.z]
	// when sh = 0
	let harmonic = [att.f_dc_0.array[i], att.f_dc_1.array[i], att.f_dc_2.array[i]]

	// when sh = 2
	// for (let j = 0; j <= 23; j++) {
	//     harmonic.push(att[`f_rest_${j}`].array[i]);
	// }
	let opacity = att.opacity.array[i]
	let scale = [att.scale_0.array[i], att.scale_1.array[i], att.scale_2.array[i]]
	let rotation = [att.rot_0.array[i], att.rot_1.array[i], att.rot_2.array[i], att.rot_3.array[i]]
	return { position, harmonic, opacity, scale, rotation }
}

async function loadOctreeGeometry(rootNode) {

	const start = performance.now()
	const sigmoid = (m1) => 1 / (1 + Math.exp(-m1))
	// Scene bouding box
	const queue = [rootNode]
	let vertexCount = 0
	let campos = [viewMatrix[2], viewMatrix[6], viewMatrix[10]];

	// first for loop to get the count of gaussians
	while (queue.length > 0) {
		const currentNode = queue.shift();
		const currentCount = currentNode.numPoints
		vertexCount += currentCount
		// set the children of the current node
		if (currentNode.children) {
			for (let cid = 0; cid < 8; cid++) {
				const child = currentNode.children[cid];
				if (child && "geometry" in child) {
					queue.push(child);
				}
			}
		}
	}
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	const buffer = new ArrayBuffer(rowLength * vertexCount);

	// second for loop to get the gaussian data 
	let glabol_i = 0
	queue.push(rootNode)
	while (queue.length > 0) {
		const currentNode = queue.shift();
		const currentCount = currentNode.numPoints
		for (let i = 0; i < currentCount; i++) {
			const positions = new Float32Array(buffer, glabol_i * rowLength, 3);
			const scales = new Float32Array(buffer, glabol_i * rowLength + 4 * 3, 3);
			const rgbas = new Uint8ClampedArray(
				buffer,
				glabol_i * rowLength + 4 * 3 + 4 * 3,
				4,
			);
			const rots = new Uint8ClampedArray(
				buffer,
				glabol_i * rowLength + 4 * 3 + 4 * 3 + 4,
				4,
			);

			let { position, harmonic, opacity, scale, rotation } = nodeView(i, currentNode)
			// Normalize quaternion
			let length2 = 0

			for (let j = 0; j < 4; j++)
				length2 += rotation[j] * rotation[j]

			const length = Math.sqrt(length2)

			rotation = rotation.map(v => (v / length) * 128 + 128)

			rots[0] = rotation[0];
			rots[1] = rotation[1];
			rots[2] = rotation[2];
			rots[3] = rotation[3];

			// Exponentiate scale
			scale = scale.map(v => Math.exp(v))

			scales[0] = scale[0];
			scales[1] = scale[1];
			scales[2] = scale[2];

			const SH_C0 = 0.28209479177387814
			rgbas[0] = (0.5 + SH_C0 * harmonic[0]) * 255;
			rgbas[1] = (0.5 + SH_C0 * harmonic[1]) * 255;
			rgbas[2] = (0.5 + SH_C0 * harmonic[2]) * 255;

			// let color = computeColorFromSH(0, position, campos, harmonic)
			// rgbas[0] = color.x * 255;
			// rgbas[1] = color.y * 255;
			// rgbas[2] = color.z * 255;

			// Activate alpha
			rgbas[3] = sigmoid(opacity) * 255;

			positions[0] = position[0];
			positions[1] = position[1];
			positions[2] = position[2];

			glabol_i++;

		}
		// set the children of the current node
		if (currentNode.children) {
			for (let cid = 0; cid < 8; cid++) {
				const child = currentNode.children[cid];
				if (child && "geometry" in child) {
					queue.push(child);
				}
			}
		}
	}

	const loadTime = `${((performance.now() - start) / 1000).toFixed(3)}s`
	console.log(`[Loader] load ${vertexCount} gaussians in ${loadTime}.`)
	return buffer;
}
let octreeFileUrl;
async function loadScene() {
	const content = await loadPointCloud(dataSource.metaData, 'geometry', octreeFileUrl)
	octreeGeometry = content.geometry
	octreeGeometryLoader = octreeGeometry.loader
	await new Promise((resolve) => setTimeout(resolve, 500));
	// wait for the geometry to load
	console.log('octreeGeometry', octreeGeometry)
	settings.renderingMode = "Octree";

	const queue = [{ node: octreeGeometry.root, level: 0 }];
	while (queue.length > 0) {
		const { node, level } = queue.shift();
		// init only handle node with level <= 4
		if (level <= 4) {
			// Load and pre-process gaussian data from .bin file
			await octreeGeometryLoader.load(node, octreeFileUrl);

			for (let cid = 0; cid < 8; cid++) {
				const child = node.children[cid];
				if (child) {
					// push child node into depth queue
					queue.push({ node: child, level: level + 1 });
				}
			}
		}
	}

	await new Promise((resolve) => setTimeout(resolve, 50));
	// init gaussian data and send to worker
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	gaussianSplats.rowLength = 3 * 4 + 3 * 4 + 4 + 4
	gaussianSplats.baseBuffer = await loadOctreeGeometry(octreeGeometry.root)
	gaussianSplats.baseVertexCount = gaussianSplats.baseBuffer.byteLength / gaussianSplats.rowLength
	// Send gaussian data to the worker
	worker.postMessage({
		buffer: gaussianSplats.baseBuffer,
		vertexCount: gaussianSplats.baseVertexCount,
	})
}
var progressTextDom = document.getElementById("progress-text");
var progressBarDom = document.getElementById("progressBar");
const chunkSize = 1024 * 1024 * 1024;//1G


const startLoadScene = (blob) => {
	octreeFileUrl = URL.createObjectURL(blob);
	progressBarDom.style.width = `100%`;
	loadScene()
}
async function getFileData() {
	const localBlob = await localforage.getItem(`${dataSource.name}_0`);
	if (localBlob) {
		const mergedBlob = await mergeChunks();
		startLoadScene(mergedBlob)
	}
	else {
		let data = await fetchSetFile(dataSource.octreeUrl)
		const mergedBlob = await mergeChunks(data);
		startLoadScene(mergedBlob)
	}
}
async function fetchSetFile(url) {
	const totalSize = await getTotalSize(url);
	const chunkCount = Math.ceil(totalSize / chunkSize);
	let chunks = []

	for (let i = 0; i < chunkCount; i++) {
		const start = i * chunkSize;
		const end = Math.min(totalSize, (i + 1) * chunkSize) - 1;

		const xhr = new XMLHttpRequest();
		xhr.open('GET', url);
		xhr.setRequestHeader('Range', `bytes=${start}-${end}`);
		xhr.responseType = 'arraybuffer';

		const promise = new Promise((resolve, reject) => {
			xhr.onload = async () => {
				if (xhr.status === 206) {
					let blob = new Blob([xhr.response]);
					let url = URL.createObjectURL(blob);
					// await storeChunkInIndexedDB(blob, i);
					resolve(blob);
				} else {
					reject(new Error(`Failed to fetch chunk: ${xhr.statusText}`));
				}
			};
			xhr.onerror = () => {
				reject(new Error('Failed to fetch chunk: Network error'));
			};
		});
		// 
		xhr.addEventListener('progress', function (event) {
			if (event.lengthComputable) {
				const progress = Math.round((event.loaded / event.total) * 100);
				progressBarDom.style.width = `${progress}%`;
				progressTextDom.innerHTML = `Load Resources ${(dataSource.size * progress / 100).toFixed(2)}Mb`;
				if (progress === 100) {
					progressTextDom.innerText = `File Loaded Successfully, Parsing...`;
				}
			} else {
				console.log('Download progress: Not computable');
			}
		});

		xhr.send();
		chunks.push(promise);
	}

	return Promise.all(chunks);
}

async function getTotalSize(url) {
	const xhr = new XMLHttpRequest();
	xhr.open('HEAD', url);
	xhr.send();
	await new Promise(resolve => xhr.onload = resolve);
	return parseInt(xhr.getResponseHeader('Content-Length'));
}

async function storeChunkInIndexedDB(chunk, index) {
	const dbName = `${dataSource.name}_${index}`;
	await localforage.setItem(dbName, chunk);
}

// 获取分片数据并合并
async function mergeChunks(data) {
	const chunks = [];
	let index = 0;
	while (true) {
		let chunk
		if (data) {
			chunk = data[index]
		}
		else {
			chunk = await localforage.getItem(`${dataSource.name}_${index}`);
		}
		if (!chunk) break;
		chunks.push(chunk);
		index++;
	}
	const blob = new Blob(chunks);
	return blob;
}

getFileData()
// fetchSetFile(dataSource.octreeUrl)

async function updateGaussianByView(viewMatrix, projectionMatrix, maxLevel, maxCount) {
	const start = performance.now()
	const queue = [{ node: octreeGeometry.root, level: 0 }];
	let ZDepthMax = 0;
	let extLoadedCount = 0;
	stopReading = false;
	gaussianSplats.extraVertexCount = 0;
	// first loop to get far dist in current view and count all extra gaussians
	while (queue.length > 0) {
		const { node, level } = queue.shift();
		if (level <= maxLevel) {
			let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, node.boundingSphere.center)
			ZDepthMax = ZDepthMax > ZDepth ? ZDepthMax : ZDepth;
			if (visibility && level > 4) extLoadedCount += node.numPoints;
		}

		// set the children of the current node
		for (let cid = 0; cid < 8; cid++) {
			const child = node.children[cid];
			if (child) {
				// push child node into depth queue
				queue.push({ node: child, level: level + 1 });
			}
		}

	}

	// second loop to get the valid count gaussians
	queue.push({ node: octreeGeometry.root, level: 0 });
	while (queue.length > 0) {
		const { node, level } = queue.shift();
		node.reading = false;
		node.visibility = false;
		if (stopReading) return;
		if (level <= maxLevel) {
			let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, node.boundingSphere.center)

			// default actLevel def by z depth
			if (extLoadedCount > maxCount) {
				const actLevel = maxLevel / Math.exp(settings.scalingBeta * ZDepth / ZDepthMax);
				node.visibility = level <= actLevel ? visibility : false;
			} else {
				node.visibility = visibility
			}

			// default level <= 4 is visible, and already preloaded
			if (level <= 4) {
				node.visibility = true;
				node.reading = true;
			}

			// Load and pre-process gaussian data from .bin file for visible node
			if (node.visibility && !node.reading) {
				await octreeGeometryLoader.load(node, octreeFileUrl);
				gaussianSplats.extraVertexCount += node.numPoints;
				if (gaussianSplats.extraVertexCount > maxCount) {
					break;
				}
			}

			// set the children of the current node
			for (let cid = 0; cid < 8; cid++) {
				const child = node.children[cid];
				if (child) {
					// push child node into depth queue
					queue.push({ node: child, level: level + 1 });
				}
			}
		}
	}

	gaussianSplats.extraBuffer = new ArrayBuffer(gaussianSplats.rowLength * gaussianSplats.extraVertexCount);
	gaussianSplats.loadedCount = 0;
	gaussianSplats.lastloadedCount = 0;

	let maxLoop = 32;
	let campos = [viewMatrix[2], viewMatrix[6], viewMatrix[10]];
	// load gaussian data into buffer
	while (true) {
		if (gaussianSplats.loadedCount >= gaussianSplats.extraVertexCount || stopReading || maxLoop-- <= 0) {
			console.log("Stop Reading Guassian Geometry!", gaussianSplats.loadedCount, gaussianSplats.extraVertexCount, maxLoop)
			break;
		}
		await readGaussianFromNode(octreeGeometry.root, gaussianSplats, campos, 0);
		// Send gaussianSplats data to the worker
		if (gaussianSplats.loadedCount > gaussianSplats.lastloadedCount + 200000 || gaussianSplats.loadedCount == gaussianSplats.extraVertexCount) {
			if (gaussianSplats.loadedCount == gaussianSplats.extraVertexCount) stopReading = true;
			let buffer = new ArrayBuffer(gaussianSplats.baseBuffer.byteLength + gaussianSplats.extraBuffer.byteLength);
			let bufferView = new Uint8Array(buffer);
			bufferView.set(new Uint8Array(gaussianSplats.baseBuffer), 0);
			bufferView.set(new Uint8Array(gaussianSplats.extraBuffer), gaussianSplats.baseBuffer.byteLength);
			worker.postMessage({
				buffer: buffer,
				vertexCount: gaussianSplats.baseVertexCount + gaussianSplats.extraVertexCount,
			})
			gaussianSplats.lastloadedCount = gaussianSplats.loadedCount;
		}
	}

	const loadTime = `${((performance.now() - start) / 1000).toFixed(3)}s`
	progressTextDom.innerHTML = ``;
	console.log(`[Loader] load ${gaussianSplats.extraVertexCount} gaussians in ${loadTime}.`)
}
// read gaussian data from octree node
async function readGaussianFromNode(node, gaussianSplats, campos, level) {
	if (node.visibility && node.loaded && !node.reading) {
		//set reading flag
		node.reading = true;
		const currentCount = node.numPoints
		for (let i = 0; i < currentCount; i++) {
			// check buffer size
			if (!isBufferLargeEnough(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3, 4)) return;
			// read into buffer
			const positions = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3);
			const scales = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3, 3);
			const rgbas = new Uint8ClampedArray(
				gaussianSplats.extraBuffer,
				gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3,
				4,
			);
			const rots = new Uint8ClampedArray(
				gaussianSplats.extraBuffer,
				gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3 + 4,
				4,
			);

			let { position, harmonic, opacity, scale, rotation } = nodeView(i, node)
			// Normalize quaternion
			let length2 = 0

			for (let j = 0; j < 4; j++)
				length2 += rotation[j] * rotation[j]

			const length = Math.sqrt(length2)

			rotation = rotation.map(v => (v / length) * 128 + 128)

			rots[0] = rotation[0];
			rots[1] = rotation[1];
			rots[2] = rotation[2];
			rots[3] = rotation[3];

			// Exponentiate scale
			scale = scale.map(v => Math.exp(v))

			scales[0] = scale[0];
			scales[1] = scale[1];
			scales[2] = scale[2];

			const SH_C0 = 0.28209479177387814
			rgbas[0] = (0.5 + SH_C0 * harmonic[0]) * 255;
			rgbas[1] = (0.5 + SH_C0 * harmonic[1]) * 255;
			rgbas[2] = (0.5 + SH_C0 * harmonic[2]) * 255;

			// let color = computeColorFromSH(0, position, campos, harmonic)
			// rgbas[0] = color.x * 255;
			// rgbas[1] = color.y * 255;
			// rgbas[2] = color.z * 255;

			// Activate alpha
			const sigmoid = (m1) => 1 / (1 + Math.exp(-m1))
			rgbas[3] = sigmoid(opacity) * 255;

			positions[0] = position[0];
			positions[1] = position[1];
			positions[2] = position[2];

			gaussianSplats.loadedCount++;
		}
	}
	// set the children of the current node
	for (let cid = 0; cid < 8; cid++) {
		const child = node.children[cid];
		if (child) {
			await readGaussianFromNode(child, gaussianSplats, level + 1)
		}
	}
}


function markCubeVisibility(viewMatrix, projectionMatrix, centerPos) {
	const x = centerPos.x;
	const y = centerPos.y;
	const z = centerPos.z;

	// Convert to camara space
	const ZDepth = x * viewMatrix[2] + y * viewMatrix[6] + z * viewMatrix[10] + viewMatrix[14];

	//Get View Projection Matrix
	const viewProjMatrix = multiply4(projectionMatrix, viewMatrix);

	// Convert to clip space
	let clipX = x * viewProjMatrix[0] + y * viewProjMatrix[4] + z * viewProjMatrix[8] + viewProjMatrix[12];
	let clipY = x * viewProjMatrix[1] + y * viewProjMatrix[5] + z * viewProjMatrix[9] + viewProjMatrix[13];
	let clipZ = x * viewProjMatrix[2] + y * viewProjMatrix[6] + z * viewProjMatrix[10] + viewProjMatrix[14];
	let clipW = x * viewProjMatrix[3] + y * viewProjMatrix[7] + z * viewProjMatrix[11] + viewProjMatrix[15];

	clipX = clipX / (clipW + 0.000001);
	clipY = clipY / (clipW + 0.000001);
	clipZ = clipZ / (clipW + 0.000001);

	let visibility = false;
	// Check if the point is within the camera's field of view
	if (ZDepth > -0.5 && clipX > -1.3 && clipX < 1.3 && clipY > -1.3 && clipY < 1.3) { //  
		visibility = true;
	}
	return { ZDepth, visibility };
}

function isBufferLargeEnough(buffer, offset, length, byteSize) {
	let requiredBytes = length * byteSize;
	return buffer.byteLength >= (offset + requiredBytes);
}

function computeColorFromSH(deg, position, campos, harmonic) {
	const SH_C0 = 0.28209479177387814;
	const SH_C1 = 0.4886025119029199;
	const SH_C2 = [1.0925484305920792,
		-1.0925484305920792,
		0.31539156525252005,
		-1.0925484305920792,
		0.5462742152960396];

	const SH_C3 = [
		-0.5900435899266435,
		2.890611442640554,
		-0.4570457994644658,
		0.3731763325901154,
		-0.4570457994644658,
		1.445305721320277,
		-0.5900435899266435
	];

	// Helper function to normalize a vector
	function normalize(vec) {
		let len = Math.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
		return { x: vec.x / len, y: vec.y / len, z: vec.z / len };
	}

	// Helper function to compute the maximum of a vector and 0
	function maxVec(vec) {
		return { x: Math.max(vec.x, 0), y: Math.max(vec.y, 0), z: Math.max(vec.z, 0) };
	}

	let pos = { x: position[0], y: position[1], z: position[2] };;
	let dir = { x: pos.x - campos[0], y: pos.y - campos[1], z: pos.z - campos[2] };
	dir = normalize(dir);

	let sh = harmonic;
	let result = { x: SH_C0 * sh[0], y: SH_C0 * sh[1], z: SH_C0 * sh[2] };

	if (deg > 0) {
		let x = dir.x;
		let y = dir.y;
		let z = dir.z;
		result = {
			x: result.x - SH_C1 * y * sh[3] + SH_C1 * z * sh[6] - SH_C1 * x * sh[9],
			y: result.y - SH_C1 * y * sh[4] + SH_C1 * z * sh[7] - SH_C1 * x * sh[10],
			z: result.z - SH_C1 * y * sh[5] + SH_C1 * z * sh[8] - SH_C1 * x * sh[11]
		};

		if (deg > 1) {
			let xx = x * x, yy = y * y, zz = z * z;
			let xy = x * y, yz = y * z, xz = x * z;
			result = {
				x: result.x + SH_C2[0] * xy * sh[12] + SH_C2[1] * yz * sh[15].x + SH_C2[2] * (2.0 * zz - xx - yy) * sh[18].x + SH_C2[3] * xz * sh[21].x + SH_C2[4] * (xx - yy) * sh[24].x,
				y: result.y + SH_C2[0] * xy * sh[13] + SH_C2[1] * yz * sh[16].y + SH_C2[2] * (2.0 * zz - xx - yy) * sh[19].y + SH_C2[3] * xz * sh[22].y + SH_C2[4] * (xx - yy) * sh[25].y,
				z: result.z + SH_C2[0] * xy * sh[14] + SH_C2[1] * yz * sh[17].z + SH_C2[2] * (2.0 * zz - xx - yy) * sh[20].z + SH_C2[3] * xz * sh[23].z + SH_C2[4] * (xx - yy) * sh[26].z
			};

			if (deg > 2) {
				// Add the third-degree SH terms
				// Similar to the second-degree terms, add the necessary calculations
			}
		}
	}
	result.x += 0.5;
	result.y += 0.5;
	result.z += 0.5;

	return maxVec(result);
}
